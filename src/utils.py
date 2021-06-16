# PYTHON PACKAGES
import matplotlib.pyplot as plt
import numpy as np
import scipy.signal as signal

# MY PACKAGES
import physical_constant_units as phy


def gaussian(x, mu, sig):
    """Returns a gaussian function"""
    return (1 / sig / np.sqrt(2 * np.pi)) * np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))


def sigmoid(x, w, tau):
    """Returns a sigmoid function"""
    return 1 / (1 + np.exp(-w * (x - tau)))


def ReLU(x, x0, m):
    """Returns the ReLU function around x0 with slope m"""
    return m * (x - x0) * (x - x0 > 0)


def flat_sine(x, x0, f):
    """Returns a flat start to a sine wave with frequency f"""
    return np.sin(2 * np.pi * f * x) * (2 * np.pi * f * x >= 2 * np.pi)


def smooth(y, box_pts=10):
    """Smoothing a signal"""
    box = np.ones(box_pts) / box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth


def l2(x):
    """Compute the l2 norm of x: ||x||_2"""
    return np.sqrt(np.sum(np.power(x, 2)))


def butter_low_pass_filter(x, order, cutoff_hz, sampling_frequency_hz, plot_filter=False):
    """Returns the filtered signal using a Butterworth low pass filter"""
    b, a = signal.butter(order, cutoff_hz / (sampling_frequency_hz / 2), btype='low', analog=False)
    filtered_x = signal.filtfilt(b, a, x)
    if plot_filter:
        # plot the filter frequency response
        w, h = signal.freqz(b, a, fs=sampling_frequency_hz)
        plt.semilogx(w, 20 * np.log10(abs(h)))
        plt.title('Butterworth filter frequency response')
        plt.xlabel('Frequency [rad / sec]')
        plt.ylabel('Amplitude [dB]')
        plt.grid(which='both', axis='both')
        plt.axvline(cutoff_hz, color='green') # cutoff frequency
        plt.show()
    return filtered_x


def psd_compare(signals_list, sampling_frequency_hz, noise_amplitude=None, names=None):
    """ Plot the Power Spectral Density (PSD) of all signals in the signals list
        in units of \sqrt(power) / \sqrt(Hz)
    """
    frequencies = []
    psd = []
    for x in signals_list:
        f, Pxx = signal.welch(x, fs=sampling_frequency_hz)
        frequencies.append(f)
        psd.append(np.sqrt(Pxx))

    fig = plt.figure(figsize=(20, 4))
    ax = plt.subplot(111)
    ax.set_title('The Power Spectral Density (PSD) plot')
    for i in range(len(psd)):
        if names is not None and len(names) == len(psd):
            ax.semilogy(frequencies[i], psd[i], label=names[i])
        else:
            ax.semilogy(frequencies[i], psd[i], label=str(i))
    if noise_amplitude is not None:
        ax.hlines(noise_amplitude, frequencies[i].min(), frequencies[i].max(), label='noise amplitude')
    ax.set_ylabel(r'PSD $\left[\sqrt{\frac{a.u.}{Hz}}\right]$')
    ax.set_xlabel('Frequency [Hz]')
    ax.legend()
    ax.grid(True)
    plt.show()

