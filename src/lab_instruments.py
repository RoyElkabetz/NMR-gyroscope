import numpy as np
import scipy as sp
import scipy.signal as sps
import matplotlib.pyplot as plt


class LIA:
    def __init__(self, lpf_params: dict, alpha: np.float = 0.,
                 amp: np.float = 1.):
        self.lpf_params = lpf_params
        self.alpha = alpha
        self.amp = amp

        """set up low pass filter parameters using a Butterworth low pass filter"""
        b, a = sps.butter(self.lpf_params['order'],
                          self.lpf_params['cutoff_hz'] / (self.lpf_params['sampling_frequency_hz'] / 2),
                          btype='low', analog=False, output='ba')
        if self.lpf_params['plot_filter']:
            # plot the filter frequency response
            w, h = sps.freqz(b, a, fs=self.lpf_params['sampling_frequency_hz'])
            plt.semilogx(w, 20 * np.log10(abs(h)))
            plt.title('Butterworth filter frequency response')
            plt.xlabel('Frequency [rad / sec]')
            plt.ylabel('Amplitude [dB]')
            plt.grid(which='both', axis='both')
            plt.axvline(self.lpf_params['cutoff_hz'], color='green')  # cutoff frequency
            plt.show()

        self.lpf_a = a
        self.lpf_b = b

    def filter_signal(self, input_signal):
        """
        Low pass filter implementation
        :param input_signal: The signal that we want to filter
        :return: the filtered signal
        """
        filtered_x = sps.filtfilt(self.lpf_b, self.lpf_a, input_signal)
        return filtered_x

    def use(self, input_signal, input_time, ref_frequency):
        """

        :param input_signal: input signal array to the lock-in
        :param input_time: input time array to the lock-in
        :param ref_frequency: the reference frequency which is used for generating the reference signal
        :return:
        """
        assert len(input_signal) == len(input_time)

        # generating the X, Y reference signals
        X_ref_signal = np.cos(2 * np.pi * ref_frequency * input_time + self.alpha)
        Y_ref_signal = np.cos(2 * np.pi * ref_frequency * input_time + self.alpha + np.pi / 2)

        # get the input signal amplitude
        amplitude_input = np.max(input_signal - np.mean(input_signal))

        # modulating the input signal
        X_modulated_signal = np.multiply(input_signal, X_ref_signal) * 2 / amplitude_input
        Y_modulated_signal = np.multiply(input_signal, Y_ref_signal) * 2 / amplitude_input

        # filter and amplify modulated signals
        X_lia = self.filter_signal(X_modulated_signal) * self.amp
        Y_lia = self.filter_signal(Y_modulated_signal) * self.amp
        R_lia = np.sqrt(X_lia ** 2 + Y_lia ** 2) * self.amp

        # set zeros in signal to numerical zero
        zeros_idx = np.abs(X_lia) == 0.
        X_lia = X_lia + zeros_idx * 1e-15

        # compute the phase difference between the input signal and reference
        Theta_lia = np.arctan(Y_lia / X_lia)

        return X_lia, Y_lia, R_lia, Theta_lia

    def scan_alpha(self):
        pass


if __name__ == "__main__":
    freq = 2
    fs = 1000
    tf = 5
    phi = np.pi / 5
    amp = 4
    t = np.linspace(0, tf, tf * fs)
    x_ref = amp * np.cos(2 * np.pi * freq * t)
    # x = amp * np.cos(2 * np.pi * freq * t + phi) * (1 + 0.1 * np.random.normal(size=len(x_ref)))
    x = amp * np.cos(2 * np.pi * freq * t + phi)
    # lpf_params = {'order': 3, 'cutoff_hz': freq * (1 / fs), 'sampling_frequency_hz': fs, 'plot_filter': True}
    lpf_params = {'order': 3, 'cutoff_hz': freq * 1e-1, 'sampling_frequency_hz': fs, 'plot_filter': False}
    my_lia = LIA(lpf_params)
    X_lia, Y_lia, R_lia, Theta_lia = my_lia.use(x, t, freq)

    plt.figure(figsize=(20, 4))
    # plt.title(f'computed phase diff: {lia_phi} [rad], vs injected: {phi}\ncomputed amplitude: {lia_amp}, vs injected {amp}')
    plt.plot(t, x, label='signal')
    plt.plot(t, X_lia, label='X_lia')
    plt.plot(t, Y_lia, label='Y_lia')
    plt.plot(t, R_lia, label='R_lia')
    plt.plot(t, Theta_lia, label='Theta_lia')
    plt.legend()
    plt.show()
