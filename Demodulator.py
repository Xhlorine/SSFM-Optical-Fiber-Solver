import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sci

class Demodulator:
    def __init__(self, bps=1e6, fs=1e8, filter=[]):
        self.setBps(bps).setFs(fs).setFilter(filter)

    # Set bps
    def setBps(self, bps):
        self.bps = bps
        return self

    # Set fs
    def setFs(self, fs):
        self.fs = fs
        return self
    
    # Set filter
    def setFilter(self, filter):
        try:
            self.Ns = int(self.fs / self.bps)
            self.filter = sci.resample(filter, self.Ns)
        except AttributeError:
            print("Please set bps and fs first.")
        return self
    
    def setSignal(self, signal):
        self.signal = signal
        return self
    
    def setCorrectBits(self, correctBits):
        self.correct = correctBits
        return self
    
    # Demodulate - signal coherence
    def demodulate(self):
        sig = np.fft.fftshift(np.fft.ifft(np.fft.fftshift(self.filter)))
        self.N = self.signal.size
        self.Ns = int(self.fs / self.bps)
        self.n = self.N // self.Ns

        self.samples = np.zeros(self.n, dtype=complex)
        for i in range(self.n):
            self.samples[i] = np.sum(self.signal[(i*self.Ns):(i+1)*self.Ns] * sig)

        self.gate = np.sum(np.abs(sig)) / 2
        self.bits = (np.abs(self.samples) > self.gate).astype(int)


    def plotConstellation(self):
        ones = self.correct == 1
        zeros = self.correct == 0
        plt.figure(num='Star plot')
        plt.subplot().set_aspect('equal')
        plt.scatter(np.real(self.samples[ones]), np.imag(self.samples[ones]), marker='o', c='r', label='1')
        plt.scatter(np.real(self.samples[zeros]), np.imag(self.samples[zeros]), marker='o', c='b', label='0')
        plt.legend()
        # plt.axvline(x=self.gate, c='g')
        plt.xlabel('Real')
        plt.ylabel('Imag')
        plt.title('Star plot')
        plt.grid(True, 'major')
        m = np.max(np.abs(plt.axis())) * 1.2
        plt.axis(np.array([-1, 1, -1,  1]) * m)
        plt.axvline(0, c='g')
        plt.axhline(0, c='g')
        plt.show()

