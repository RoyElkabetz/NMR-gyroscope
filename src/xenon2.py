# PYTHON PACKAGES
from scipy.integrate import odeint
from scipy.linalg import inv
import matplotlib.pyplot as plt
import numpy as np

# MY PACKAGES
import physical_constant_units as phy
import utils


class Xenon:
    def __init__(self, gamma, t1, t2, K0=np.array([0.05, 0.03, 0.95]), ts=500, dt=1, name='129'):

        self.name = name
        # Physical parameters
        self.gamma = gamma      # The xenon gyromagnetic ration   [rad  s^-1 T^-1]

        # Decays
        self.t1 = t1            # T1 of xenon                     [s]
        self.t2 = t2            # T2 of xenon                     [s]
        self.gamma1 = 1. / t1   #                                 [s^-1]
        self.gamma2 = 1. / t2   #                                 [s^-1]

        # solver parameters
        self.ts = ts                                                # solver time frame [s]
        self.dt = dt                                                # solver time steps [s]
        self.t_steps = int(ts // dt)                                # solver number of steps
        self.time_vec = np.linspace(0, self.ts, self.t_steps)       # time vector [s]
        self.solver_done = False                                    # boolean for state of solver

        # Spin polarization
        self.Kt = np.zeros((self.t_steps, 3))
        self.Kt[0, :] = K0
        self.Ks = np.zeros((self.t_steps, 3))
        self.Rse = np.array([0, 0, 0])
        self.Kt_perp = None
        self.Ks_perp = None
        self.phase_perp = None

        # Bloch matrix
        self.M = np.zeros((3, 3))

        # boolean params
        self.drive = True

    def set_spin_exchange_amp(self, rse):
        """Set spin exchange pumping"""
        self.Rse = rse

    def set_bloch_matrix(self, environment):
        """Constructing the Bloch matrix of the dynamics"""
        i = environment.i
        M = np.array(
            [
              [-self.gamma2 , 0              , 0             ],
              [0            , -self.gamma2   , 0             ],
              [0            , 0              , -self.gamma1  ]
              ]
        )
        M12 = self.gamma * (environment.B0[i] + environment.Bnoise[i]) + environment.wr[i]
        if self.drive:
            M[0, 2] = -environment.Ad_y[i] / 2.
            M[1, 2] = environment.Ad_x[i] / 2.
            M[2, 0] = environment.Ad_y[i] / 2.
            M[2, 1] = -environment.Ad_x[i] / 2.
            M12 = M12 - environment.wd_x[i] - environment.wd_y[i]

        M[0, 1] = M12
        M[1, 0] = -M12
        self.M = M

    def solve_steady_state(self, i):
        self.Ks[i] = -inv(self.M) @ self.Rse

    def bloch_equations(self, K, t):
        """Bloch dynamics model"""
        dK_dt = np.matmul(self.M, K) + self.Rse
        return dK_dt

    def update_solver_time_frame(self, ts, dt):
        """Updating the solver parameters"""
        # solver parameters
        self.ts = ts                                            # solver time frame [s]
        self.dt = dt                                            # solver time steps [s]
        self.t_steps = int(ts // dt)                            # solver number of steps
        self.time_vec = np.linspace(0, self.ts, self.t_steps)   # time vector [s]

    def solve_dynamics(self, environment):
        """Solving the Bloch equations"""
        ts_frame = np.linspace(0, self.dt, 2)  # single time frame for solver
        for i in range(1, self.t_steps):
            environment.set_step(i - 1)
            self.set_bloch_matrix(environment=environment)
            self.solve_steady_state(i - 1)
            Kt_temp = self.Kt[i - 1, :]
            Kt_temp = odeint(self.bloch_equations, Kt_temp, ts_frame)
            self.Kt[i, :] = Kt_temp[-1, :]
        environment.set_step(self.t_steps - 1)
        self.solve_steady_state(self.t_steps - 1)
        self.solver_done = True

    def compute_perpendicular_values(self):
        """Compute perpendicular polarization magnitude and phase with respect to the drive"""
        self.Kt_perp = np.sqrt(self.Kt[:, 0] ** 2 + self.Kt[:, 1] ** 2)
        self.Ks_perp = np.sqrt(self.Ks[:, 0] ** 2 + self.Ks[:, 1] ** 2)
        self.phase_perp = np.arctan(self.Kt[:, 1] / self.Kt[:, 0])

    def display_params(self):
        print('===================================================================')
        print(f'| Xenon {self.name}:')
        print(f'| ----------')
        print(f'| gyromagnetic ratio:     {self.gamma}')
        print(f'| T1:                     {self.t1}')
        print(f'| T2:                     {self.t2}')
        print(f'| K0:                     {self.Kt[0, :]}')
        print(f'| Kt:                     {self.Kt[-1, :]}')
        print(f'| K steady:               {self.Ks[-1, :]}')
        print('===================================================================')

    def plot_results(self, environment, ti=0):
        assert self.solver_done
        self.compute_perpendicular_values()
        ts = self.time_vec
        # Plot the spin solution
        plt.rcParams.update({'font.size': 12})  # increase the font size
        plt.rcParams['axes.facecolor'] = 'white'
        plt.rcParams["figure.facecolor"] = 'lightyellow'
        plt.rcParams["figure.facecolor"] = 'lightyellow'
        plt.rcParams["lines.linewidth"] = 2.5

        fig = plt.figure(figsize=(20, 14))
        ax1 = plt.subplot(3, 2, (1, 2))
        ax1.set_title('Kx, Ky, Kz polarizations')
        ax1.set_xlabel("time [s]")
        ax1.set_ylabel("$K$")
        ax1.plot(ts, self.Kt[:, 0], label='$K_x$', color='tab:blue')
        ax1.plot(ts, self.Kt[:, 1], label='$K_y$', color='tab:orange')
        ax1.plot(ts, self.Kt[:, 2], label='$K_z$', color='tab:green')
        ax1.plot(ts, self.Ks[:, 0], '--', label='$K_x$ - steady state solution', color='tab:blue')
        ax1.plot(ts, self.Ks[:, 1], '--', label='$K_y$ - steady state solution', color='tab:orange')
        ax1.plot(ts, self.Ks[:, 2], '--', label='$K_z$ - steady state solution', color='tab:green')
        ax1.legend()
        ax1.grid(True)

        ax2 = plt.subplot(3, 2, 3)
        ax2.plot(ts, self.Kt_perp, label='$|K_{\perp}|$', color='tab:blue')
        ax2.plot(ts, self.Ks_perp, '--', label='$|K_{\perp}|$ steady state solution', color='tab:blue')
        ax2.set_xlabel("time [s]")
        ax2.set_ylabel("$K$")
        ax2.grid(True)
        ax2.legend()
        ax2.set_title('Perpendicular polarization')

        ax3 = plt.subplot(3, 2, 4)
        ax3.plot(ts, self.phase_perp * phy.R2D, label='$\phi = arctan(K_y/K_x)$')
        ax3.set_xlabel("time [s]")
        ax3.set_ylabel("$\phi$ [degree]")
        ax3.grid(True)
        ax3.legend()
        ax3.set_title('Phase difference with respect to drive')

        world_rotation = -self.phase_perp * self.gamma2 + self.gamma * environment.B0 - environment.wd_y
        ax4 = plt.subplot(3, 2, (5, 6))
        ax4.plot(ts[ts > ti], world_rotation[ts > ti], label='$\Omega_r$ - calculated', color='tab:green')
        ax4.plot(ts[ts > ti], environment.wr[ts > ti], '--', label='$\Omega_r$ - true', color='tab:green')
        ax4.set_xlabel("time [s]")
        ax4.set_ylabel("$\omega$ [rad/s]")
        ax4.grid(True)
        ax4.legend(loc='upper right')
        ax4.set_title('World rotation')
        ax5 = ax4.twinx()
        ax5.plot(ts[ts > ti], world_rotation[ts > ti] - environment.wr[ts > ti],
                 label='$\Omega_r^{calc} - \Omega_r^{true}$', color='red')
        ax5.legend(loc='lower right')
        ax5.set_ylabel('error', color='red')
        ax5.tick_params(axis='y', labelcolor='red')

        plt.tight_layout
        plt.show()


