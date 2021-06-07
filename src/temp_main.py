# PYTHON PACKAGES
from scipy.integrate import odeint
from scipy.linalg import inv
import matplotlib.pyplot as plt
import numpy as np

# MY PACKAGES
import physical_constant_units as phy
import environment2 as env
import xenon2 as xe
import utils

# solver parameters
t_final = 500
dt = 1
steps = int(t_final // dt)
ts1 = np.linspace(0, dt, 2)
ts = np.linspace(0, t_final, steps)

# world rotation parameters
#wr = 0.01 * utils.sigmoid(ts, 1, 100)                          # rad / s
wr = 10 * utils.gaussian(ts, 200, 50)                          # rad / s

# Xenon parameters
g129 = phy.G129                        # gyromagnetic ratio [rad  s^-1  T^-1]
t1 = 30                                # s
t2 = 8                                 # s
Rse = np.array([0, 0, 0.1]) * t1       # |K| / s

# Environment parameters
B0 = 1 * phy.G2T * np.ones_like(ts)                             # Tesla
Bnoise = 1e-6 * phy.G2T * utils.smooth(np.random.randn(len(ts)))   # Tesla
Ad_y = 2 * np.sqrt((1 / t1) * (1 / t2)) * np.ones_like(ts)      # rad / s
wd_y = g129 * 1 * phy.G2T * np.ones_like(ts)                    # rad / s
Ad_x = np.zeros_like(ts)                                        # rad / s
wd_x = np.zeros_like(ts)                                        # rad / s

# initialize Xenon
Xe129 = xe.Xenon(gamma=g129, t1=t1, t2=t2, K0=np.array([0.0259, 0.02, 0.3]), ts=t_final, dt=dt, name='129')
Xe129.set_spin_exchange_amp(Rse)
Xe129.display_params()

# initialize Environment
env129 = env.Environment()
env129.set_state(wr=wr, B0=B0, Bnoise=Bnoise, Ad_y=Ad_y, wd_y=wd_y, Ad_x=Ad_x, wd_x=wd_x)
env129.display_params()

# plot world rotation
fig = plt.figure(figsize=(12, 3))
ax = plt.subplot()
ax.plot(ts, Xe129.gamma * B0, label='$\omega_0$')
ax.plot(ts, Xe129.gamma * (B0 + Bnoise), label='$\omega_0 + \omega_{noise}$')
ax.set_ylabel('$\hat{z}$ magnetic rotation [rad / s]')
ax.set_xlabel('time [s]')
ax.legend()
ax1 = ax.twinx()
ax1.plot(ts, wr, label='$\Omega_r$', color='green')
ax1.set_ylabel('$\hat{z}$ world rotation [rad / s]')
ax1.legend()
plt.show()


# run solver and save dynamics
Xe129.solve_dynamics(env129)
Xe129.plot_results(env129)
