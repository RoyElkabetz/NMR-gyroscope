# PYTHON PACKAGES
from scipy.integrate import odeint
from scipy.linalg import inv
import matplotlib.pyplot as plt
import numpy as np

# MY PACKAGES
import physical_constant_units as phy
import environment as env
import xenon as xe
import utils


# x = np.linspace(0, 100, 1000)
# x0 = 30
# m = 0.1
# y = utils.ReLU(x, x0, m)
#
# plt.figure()
# plt.plot(x, y)
# plt.show()


# # solver parameters
# t_final = 500
# dt = 1
# steps = int(t_final // dt)
# ts1 = np.linspace(0, dt, 2)
# ts = np.linspace(0, t_final, steps)
#
# # world rotation parameters
# #wr = 0.01 * utils.sigmoid(ts, 1, 100)                          # rad / s
# wr = 10 * utils.gaussian(ts, 200, 50)                          # rad / s
#
# # Xenon parameters
# g129 = phy.G129                        # gyromagnetic ratio [rad  s^-1  T^-1]
# t1 = 30                                # s
# t2 = 8                                 # s
# Rse = np.array([0, 0, 0.1]) * t1       # |K| / s
#
# # Environment parameters
# B0 = 1 * phy.G2T * np.ones_like(ts)                             # Tesla
# Bnoise = 1e-6 * phy.G2T * utils.smooth(np.random.randn(len(ts)))   # Tesla
# Ad_y = 2 * np.sqrt((1 / t1) * (1 / t2)) * np.ones_like(ts)      # rad / s
# wd_y = g129 * 1 * phy.G2T * np.ones_like(ts)                    # rad / s
# Ad_x = np.zeros_like(ts)                                        # rad / s
# wd_x = np.zeros_like(ts)                                        # rad / s
#
# # initialize Xenon
# Xe129 = xe.Xenon(gamma=g129, t1=t1, t2=t2, K0=np.array([0.0259, 0.02, 0.3]), ts=t_final, dt=dt, name='129')
# Xe129.set_spin_exchange_amp(Rse)
# Xe129.display_params()
#
# # initialize Environment
# env129 = env.Environment()
# env129.set_state(wr=wr, B0=B0, Bnoise=Bnoise, Ad_y=Ad_y, wd_y=wd_y, Ad_x=Ad_x, wd_x=wd_x)
# env129.display_params()
#
# # plot world rotation
# fig = plt.figure(figsize=(12, 3))
# ax = plt.subplot()
# ax.plot(ts, Xe129.gamma * B0, label='$\omega_0$')
# ax.plot(ts, Xe129.gamma * (B0 + Bnoise), label='$\omega_0 + \omega_{noise}$')
# ax.set_ylabel('$\hat{z}$ magnetic rotation [rad / s]')
# ax.set_xlabel('time [s]')
# ax.legend()
# ax1 = ax.twinx()
# ax1.plot(ts, wr, label='$\Omega_r$', color='green')
# ax1.set_ylabel('$\hat{z}$ world rotation [rad / s]')
# ax1.legend()
# plt.show()
#
#
# # run solver and save dynamics
# Xe129.solve_dynamics(env129)
# Xe129.plot_results(env129)


# Signal parameters
frequency = 2                           # sine wave frequency [Hz]
sampling_frequency = 20 * frequency     # sampling frequency [Hz]
t_f = 30                                # final time of time vector [s]
dt = 1 / sampling_frequency             # time step for simulation [s]
N = np.int(t_f // dt)                   # sampling points in time vector
t = np.linspace(0, t_f, N)              # time vector  [s]
signal_amplitude = 10                   # signal amplitude [amplitude]

# Noise parameters
cutoff = 5                             # noise cutoff frequency [Hz]
order = 10                             # order of the low pass filter
noise_amplitude = 0.067                # noise amplitude [amplitude / sqrt(Hz)]
noise_power = np.power(noise_amplitude, 2) * sampling_frequency / 2 # noise power [power / Hz]

# The noise
noise = np.random.normal(scale=np.sqrt(noise_power), size=t.shape)
filtered_noise = utils.butter_low_pass_filter(noise,
                                              order,
                                              cutoff,
                                              sampling_frequency,
                                              plot_filter=True) # filtered noise


# The signals
signal = signal_amplitude * np.sin(2 * np.pi * frequency * t)
noisy_signal = signal + noise
filtered_noisy_signal = signal + filtered_noise
smoothed_signal = utils.smooth(noisy_signal, 10)  # smoothing the signal with mooving average


# plot all the signals
fig = plt.figure(figsize=(20, 4))
ax = plt.subplot(111)
ax.plot(t, signal, label='signal')
ax.plot(t, noisy_signal, label='noisy signal')
ax.plot(t, filtered_noisy_signal, label='filtered noisy signal')
ax.plot(t, smoothed_signal, label='smoothed signal')
ax.set_xlabel('Time [s]')
ax.set_ylabel('Amplitude [a.u.]') # arbitrary units (units of amplitude)
ax.legend()
ax.grid(True)
plt.show()


signals_list = [noise, signal, noisy_signal, filtered_noisy_signal, smoothed_signal]
names = ['noise', 'signal', 'noisy signal', 'filtered noisy signal', 'smoothed signal']

# plot PSD of all signals
utils.psd_compare(signals_list, sampling_frequency, noise_amplitude=noise_amplitude, names=names)



