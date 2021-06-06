# PYTHON PACKAGES
from scipy.integrate import odeint
from scipy.linalg import inv
import matplotlib.pyplot as plt
import numpy as np

# MY PACKAGES
import physical_constant_units as phy
import utils


class Xenon:
    def __init__(self, gamma, t1, t2, K0=np.array([0.05, 0.03, 0.95]), name='129'):

        self.name = name
        # Physical parameters
        self.gamma = gamma      # The xenon gyromagnetic ration   [rad  s^-1 T^-1]

        # Decays
        self.t1 = t1            # T1 of xenon                     [s]
        self.t2 = t2            # T2 of xenon                     [s]
        self.gamma1 = 1. / t1   #                                 [s^-1]
        self.gamma2 = 1. / t2   #                                 [s^-1]

        # Spin polarization
        self.K0 = K0
        self.Kt = K0
        self.Ksteady = None
        self.Rse = np.array([0, 0, 0])

        # Bloch matrix
        self.M = np.zeros((3, 3))

        # boolean params
        self.drive = True

    def set_spin_exchange_amp(self, rse):
        """Set spin exchange pumping"""
        self.Rse = rse

    def set_bloch_matrix(self, environment):
        """Constructing the Bloch matrix of the dynamics"""
        M = np.array(
            [
              [-self.gamma2 , 0              , 0             ],
              [0            , -self.gamma2   , 0             ],
              [0            , 0              , -self.gamma1  ]
              ]
        )
        M12 = self.gamma * (environment.B0 + environment.Bnoise) + environment.wr
        if self.drive:
            M[0, 2] = -environment.Ad_y / 2.
            M[1, 2] = environment.Ad_x / 2.
            M[2, 0] = environment.Ad_y / 2.
            M[2, 1] = -environment.Ad_x / 2.
            M12 = M12 - environment.wd_x - environment.wd_y

        M[0, 1] = M12
        M[1, 0] = -M12
        self.M = M

    def bloch_equations(self, K, t):
        """Bloch dynamics model"""
        dK_dt = np.matmul(self.M, K) + self.Rse
        return dK_dt

    def solve_dynamics(self, ts):
        Kt = self.K0
        Kt = odeint(self.bloch_equations, Kt, ts)
        self.Kt = Kt[-1, :]
        return Kt

    def solve_steady_state(self):
        self.Ksteady = -inv(self.M) @ self.Rse

    def display_params(self):
        print('===================================================================')
        print(f'| Xenon {self.name}:')
        print(f'| ----------')
        print(f'| gyromagnetic ratio:     {self.gamma}')
        print(f'| T1:                     {self.t1}')
        print(f'| T2:                     {self.t2}')
        print(f'| K0:                     {self.K0}')
        print(f'| Kt:                     {self.Kt}')
        print(f'| K steady:               {self.Ksteady}')
        print('===================================================================')


