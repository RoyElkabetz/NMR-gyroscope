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
