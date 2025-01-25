from typing import Literal
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sci

class Demodulator:
    def __init__(self, bps=1e6, fs=1e8, filter=[], padding=0):
        self.setBps(bps).setFs(fs).setFilter(filter).setPadding(padding)

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
    
    # Set padding
    def setPadding(self, padding):
        self.padding = padding
        return self
    
    def setSignal(self, signal):
        self.signal = signal
        return self
    
    def setCorrectBits(self, correctBits, method: Literal['OOK', 'QPSK', '16QAM', '32QAM', '64QAM']):
        self.correct = correctBits
        self.method = method
        return self
    
    # Demodulate - signal coherence
    def demodulate(self):
        sig = np.fft.fftshift(np.fft.ifft(np.fft.fftshift(self.filter)))
        self.N = self.signal.size
        self.Ns = int(self.fs / self.bps)
        self.n = self.N // self.Ns

        self.samples = np.zeros(self.n-self.padding*2, dtype=complex)
        self.bits = np.zeros(self.n-self.padding*2, dtype=int)
        for i in range(self.n-self.padding*2):
            self.samples[i] = np.sum(self.signal[((i+self.padding)*self.Ns):(i+self.padding+1)*self.Ns] * sig)
            self.bits[i] = np.argmin(np.abs(symbols[self.method] - self.samples[i]))

    def plotConstellation(self):
        plt.figure(num='Star plot')
        plt.subplot().set_aspect('equal')
        plt.axvline(0, c='g')
        plt.axhline(0, c='g')
        for i in range(symbols[self.method].size):
            mask = self.correct == i
            plt.scatter(np.real(self.samples[mask]), np.imag(self.samples[mask]), marker='o', label=str(i), s=25)
        # plt.legend()
        # plt.axvline(x=self.gate, c='g')
        plt.xlabel('Real')
        plt.ylabel('Imag')
        plt.title('Constellation')
        plt.grid(True, 'major')
        m = np.max(np.abs(np.hstack((np.real(self.samples), np.imag(self.samples))))) * 1.2
        plt.axis(np.array([-1, 1, -1,  1]) * m)
        plt.show()


# Symbol mapping
# By default, the least energy of all symbols is sqrt(2).
OOK = np.array([0, 1], dtype=complex)
_temp = np.array([-1, 1], dtype=complex)
QPSK = (_temp + 1j *  _temp[:, None]).flatten()
_temp = np.array([-3, -1, 1, 3], dtype=complex)
QAM16 = (_temp + 1j *  _temp[:, None]).flatten()
_temp = np.array([-5, -3, -1, 1, 3, 5], dtype=complex)
QAM32 = (_temp + 1j *  _temp[:, None]).flatten()
QAM32 = QAM32[(np.abs(QAM32) < 5*1.414)]
_temp = np.array([-7, -5, -3, -1, 1, 3, 5, 7], dtype=complex)
QAM64 = (_temp + 1j *  _temp[:, None]).flatten()

symbols = {'OOK': OOK, 'QPSK': QPSK, '16QAM': QAM16, '32QAM': QAM32, '64QAM': QAM64}
del OOK, _temp, QPSK, QAM16, QAM32, QAM64

def symbolMapping(bits: list | np.ndarray,
                  method: Literal['OOK', 'QPSK', '16QAM', '32QAM', '64QAM']='OOK',
                  maxEnergy: float | np.number=None,
                  leastEnergy: float | np.number=None,
                  averageEnergy: float | np.number=None):
    alphabet = symbols[method]

    if maxEnergy is not None:
        factor = np.sqrt(maxEnergy) / np.max(np.abs(alphabet))
    elif leastEnergy is not None:
        factor = np.sqrt(leastEnergy) / np.min(np.abs(alphabet))
    elif averageEnergy is not None:
        factor = np.sqrt(averageEnergy) / np.mean(np.abs(alphabet))
    else:
        factor = 1.0

    return np.array(alphabet[bits] * factor)