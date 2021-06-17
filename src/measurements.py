# PYTHON PACKAGES
from scipy.integrate import odeint
from scipy.linalg import inv
import matplotlib.pyplot as plt
import numpy as np
from tqdm.autonotebook import tqdm


# MY PACKAGES
import physical_constant_units as phy
import environment as env
import xenon as xe
import utils


def single_species_Open_Loop_bandwidth_simualtion(gyromagnetic, t1, t2, wr_amp=0.01, B0_amp=1e-6, Bnoise_amp=0, noise_cutoff_hz=0.1, filter_order=2, num_periods=2, points_in_period=1000, freq_list=None, plot_results=True, get_values=False):

    if freq_list is None:
        estimated_bandwidth = 1 / t2 / np.pi
        freq_list = np.logspace(np.log10(estimated_bandwidth) - 2, np.log10(estimated_bandwidth) + 2, 30)

    phase_diff = np.zeros_like(freq_list)
    amplitude_ratio = np.zeros_like(freq_list)

    for i, f in enumerate(tqdm(freq_list)):
        # solver parameters
        freq = f                        # [Hz]
        period = 1 / freq               # [s]
        t_final = num_periods * period  # [s]
        if t_final < 10 * t2:
            t_final = 10 * t2
        dt = period / points_in_period  # [s]
        sampling_frequency = 1 / dt     # [Hz]
        steps = int(t_final // dt)
        ts1 = np.linspace(0, dt, 2)
        ts = np.linspace(0, t_final, steps)
        Bnoise = np.zeros_like(ts)
        if Bnoise_amp != 0:
            noise = utils.get_white_noise(Bnoise_amp * phy.G2T, sampling_frequency, ts)
            Bnoise = utils.butter_low_pass_filter(noise, filter_order, noise_cutoff_hz, sampling_frequency)  # Tesla

        # world rotation
        wr = wr_amp * np.sin(2 * np.pi * freq * ts)             # rad / s

        # Xenon parameters
        Rse = np.array([0, 0, 0.1]) * t1       # |K| / s

        # Environment parameters
        B0 = B0_amp * phy.G2T * np.ones_like(ts)                             # Tesla
        Ad_y = 2 * np.sqrt((1 / t1) * (1 / t2)) * np.ones_like(ts)      # rad / s
        wd_y = gyromagnetic * B0_amp * phy.G2T * np.ones_like(ts)            # rad / s
        Ad_x = np.zeros_like(ts)                                        # rad / s
        wd_x = np.zeros_like(ts)                                        # rad / s

        # initialize Environment
        my_env = env.Environment()
        my_env.set_state(wr=wr, B0=B0, Bnoise=Bnoise, Ad_y=Ad_y, wd_y=wd_y, Ad_x=Ad_x, wd_x=wd_x)

        # initialize Xenon
        my_Xe = xe.Xenon(gamma=gyromagnetic, t1=t1, t2=t2, K0=np.array([0.0259, 0.02, 0.3]), ts=t_final, dt=dt)
        my_Xe.set_spin_exchange_amp(Rse)
        my_Xe.set_bloch_matrix(my_env)
        my_Xe.init_with_steady_state()

        # run solver and save dynamics
        my_Xe.solve_dynamics(my_env)
        my_Xe.compute_perpendicular_values()

        # computing the world rotation from xenon measurements
        world_rotation = -my_Xe.phase_perp * my_Xe.gamma2 + my_Xe.gamma * my_env.B0 - my_env.wd_y
        amplitude_ratio[i] = np.max(world_rotation[ts > 8 * t2]) / np.max(wr[ts > 8 * t2])
        phase_diff[i] = np.arccos(np.dot(wr[ts > 8 * t2], world_rotation[ts > 8 * t2]) / utils.l2(wr[ts > 8 * t2]) / utils.l2(world_rotation[ts > 8 * t2])) / np.pi * 180

    if plot_results:
        # plot results
        fig = plt.figure(figsize=(12, 8))
        ax = plt.subplot(212)
        ax.semilogx(freq_list, phase_diff, 'o', label='$\Delta\phi$')
        ax.vlines(1 / t2 / np.pi, ymin=0, ymax=90,
                  color='orange', label=r'$\frac{1}{\pi T_2}$ ')
        ax.set_xlabel('Frequency [Hz]')
        ax.set_ylabel('Phase [rad]')
        ax.grid(True)
        ax.legend()

        ax1 = plt.subplot(211)
        ax1.set_title('Single species Open-Loop bandwidth simulation')
        ax1.semilogx(freq_list, 10 * np.log10(amplitude_ratio), 'o', label=r'$\frac{|\Omega_r^{Calc}|}{|\Omega_r^{True}|}$')
        ax1.vlines(1 / t2 / np.pi, ymin=np.min(10 * np.log10(amplitude_ratio)), ymax=0,
                   color='orange', label=r'$\frac{1}{\pi T_2}$')
        ax1.hlines(-3, xmin=freq_list[0], xmax=freq_list[-1],
                   color='red', label=r'$-3 [dB]$')
        ax1.set_ylabel('Magnitude [dB]')
        ax1.set_xlabel('Frequency [Hz]')
        ax1.grid(True)
        ax1.legend()
        plt.tight_layout
        plt.show()

    if get_values:
        return freq_list, phase_diff, amplitude_ratio


def single_species_Open_Loop_dynamic_range_simulation(gyromagnetic, t1, t2, wr_amp, B0_amp=1e-6, Bnoise_amp=0, noise_cutoff_hz=0.1, filter_order=2, dt=1, t_final=1000, plot_results=True, get_values=False):
    # world rotation parameters
    wr_measurements = np.zeros_like(wr_amp)

    # solver parameters
    steps = int(t_final // dt)
    sampling_frequency = 1. / dt
    ts = np.linspace(0, t_final, steps)
    Bnoise = np.zeros_like(ts)
        
    # solver
    for i, amp in enumerate(tqdm(wr_amp)):
        # world rotation
        wr = amp * utils.sigmoid(ts, 1, 100)  # rad / s

        # Xenon parameters
        Rse = np.array([0, 0, 0.1]) * t1  # |K| / s

        # Environment parameters
        B0 = B0_amp * phy.G2T * np.ones_like(ts)  # Tesla
        if Bnoise_amp != 0:
            noise = utils.get_white_noise(Bnoise_amp * phy.G2T, sampling_frequency, ts)
            Bnoise = utils.butter_low_pass_filter(noise, filter_order, noise_cutoff_hz, sampling_frequency)
        Ad_y = 2 * np.sqrt((1 / t1) * (1 / t2)) * np.ones_like(ts)  # rad / s
        wd_y = gyromagnetic * B0_amp * phy.G2T * np.ones_like(ts)  # rad / s
        Ad_x = np.zeros_like(ts)  # rad / s
        wd_x = np.zeros_like(ts)  # rad / s

        # initialize Environment
        my_env = env.Environment()
        my_env.set_state(wr=wr, B0=B0, Bnoise=Bnoise, Ad_y=Ad_y, wd_y=wd_y, Ad_x=Ad_x, wd_x=wd_x)

        # initialize Xenon
        my_Xe = xe.Xenon(gamma=gyromagnetic, t1=t1, t2=t2, K0=np.array([0.0259, 0.02, 0.3]), ts=t_final, dt=dt)
        my_Xe.set_spin_exchange_amp(Rse)
        my_Xe.set_bloch_matrix(my_env)
        my_Xe.init_with_steady_state()

        # run solver and save dynamics
        my_Xe.solve_dynamics(my_env)
        my_Xe.compute_perpendicular_values()

        # computing the world rotation from xenon measurements
        world_rotation = -my_Xe.phase_perp * my_Xe.gamma2 + my_Xe.gamma * my_env.B0 - my_env.wd_y
        wr_measurements[i] = world_rotation[-1]

    if plot_results:
        # plot dynamic range results
        fig = plt.figure(figsize=(12, 8))
        ax = plt.subplot()
        ax.set_title('World rotation dynamic range simulation')
        ax.loglog(wr_amp, wr_measurements, 'o', label='simulation')
        ax.loglog(wr_amp, wr_amp, '--', label='perfect match')
        ax.set_ylabel('$\Omega_r^{calc}$ [rad / s]')
        ax.set_xlabel('$\Omega_r^{True}$ [rad / s]')
        ax.vlines(1 / np.sqrt(t1 * t2), ymin=np.min(wr_measurements), ymax=np.max(wr_measurements),
                  color='red', label=r'$\frac{1}{\sqrt{T_1 * T_2}}$')
        if Bnoise_amp != 0:
            ax.vlines(np.abs(gyromagnetic * Bnoise_amp * phy.G2T), ymin=np.min(wr_measurements), ymax=np.max(wr_measurements),
                      color='black', label=r'$\gamma |B^{noise}|$')
        ax.grid(True)
        ax.legend()
        plt.tight_layout
        plt.show()

    if get_values:
        return wr_amp, wr_measurements
