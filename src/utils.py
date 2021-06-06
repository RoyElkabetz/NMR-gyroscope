# PYTHON PACKAGES
import matplotlib.pyplot as plt
import numpy as np

# MY PACKAGES
import physical_constant_units as phy


def gaussian(x, mu, sig):
    """Returns a gaussian function"""
    return (1 / sig / np.sqrt(2 * np.pi)) * np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))


def sigmoid(x, w, tau):
    """Returns a sigmoid function"""
    return 1 / (1 + np.exp(-w * (x - tau)))


def smooth(y, box_pts=10):
    """Smoothing a signal"""
    box = np.ones(box_pts) / box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth


def plot_results(ts, Kt, Ks, xenon, environment, wr):
    # Plot the spin solution
    plt.rcParams.update({'font.size': 12})  # increase the font size
    plt.rcParams['axes.facecolor'] = 'white'
    plt.rcParams["figure.facecolor"] = 'lightyellow'
    plt.rcParams["figure.facecolor"] = 'lightyellow'
    plt.rcParams["lines.linewidth"] = 2.5

    fig = plt.figure(figsize=(20, 14))

    ax1 = plt.subplot(3, 2, (1, 2))
    ax1.set_title('Kx, Ky, Kz polarizations' )
    ax1.set_xlabel("time [s]")
    ax1.set_ylabel("$K$")
    ax1.plot(ts, Kt[:, 0], label='$K_x$', color='tab:blue')
    ax1.plot(ts, Kt[:, 1], label='$K_y$', color='tab:orange')
    ax1.plot(ts, Kt[:, 2], label='$K_z$', color='tab:green')
    ax1.plot(ts, Ks[:, 0], '--', label='$K_x$ - steady state solution', color='tab:blue')
    ax1.plot(ts, Ks[:, 1], '--', label='$K_y$ - steady state solution', color='tab:orange')
    ax1.plot(ts, Ks[:, 2], '--',  label='$K_z$ - steady state solution', color='tab:green')
    ax1.legend()
    ax1.grid(True)


    K_perp = np.sqrt(Kt[:, 0] ** 2 + Kt[:, 1] ** 2)
    K_perp_steady = np.sqrt(Ks[:, 0] ** 2 + Ks[:, 1] ** 2)
    ax2 = plt.subplot(3, 2, 3)
    ax2.plot(ts, K_perp, label='$|K_{\perp}|$', color='tab:blue')
    ax2.plot(ts, K_perp_steady, '--', label='$|K_{\perp}|$ steady state solution', color='tab:blue')
    ax2.set_xlabel("time [s]")
    ax2.set_ylabel("$K$")
    ax2.grid(True)
    ax2.legend()
    ax2.set_title('Perpendicular polarization')

    phi = np.arctan(Kt[:, 1] / Kt[:, 0])
    ax3 = plt.subplot(3, 2, 4)
    ax3.plot(ts, phi * phy.R2D, label='$\phi = arctan(K_y/K_x)$')
    ax3.set_xlabel("time [s]")
    ax3.set_ylabel("$\phi$ [degree]")
    ax3.grid(True)
    ax3.legend()
    ax3.set_title('Phase difference with respect to drive')

    ti = 0
    world_rotation = -phi * xenon.gamma2 + xenon.gamma * environment.B0 - environment.wd_y
    ax4 = plt.subplot(3, 2, (5, 6))
    ax4.plot(ts[ts > ti], world_rotation[ts > ti], label='$\Omega_r$ - calculated', color='tab:green')
    ax4.plot(ts[ts > ti], wr[ts > ti], '--', label='$\Omega_r$ - true', color='tab:green')
    ax4.set_xlabel("time [s]")
    ax4.set_ylabel("$\omega$ [rad/s]")
    ax4.grid(True)
    ax4.legend(loc='upper right')
    ax4.set_title('World rotation')
    ax5 = ax4.twinx()
    ax5.plot(ts[ts > ti], world_rotation[ts > ti] - wr[ts > ti] ,label='$\Omega_r^{calc} - \Omega_r^{true}$', color='red')
    ax5.legend(loc='lower right')
    ax5.set_ylabel('error', color='red')
    ax5.tick_params(axis='y', labelcolor='red')

    plt.tight_layout
    plt.show()
