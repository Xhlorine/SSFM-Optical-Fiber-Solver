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

        self.samples = np.zeros(shape=self.correct.shape, dtype=complex)
        self.bits = np.zeros(shape=self.correct.shape, dtype=int)
        for i in range(self.n-self.padding*2):
            self.samples[i] = np.sum(self.signal[((i+self.padding)*self.Ns):(i+self.padding+1)*self.Ns] * sig)
            self.bits[i] = np.argmin(np.abs(symbols[self.method] - self.samples[i]))
        return self

    def plotConstellation(self, legend:bool=True):
        plt.figure(num='Constellation')
        plt.subplot().set_aspect('equal')
        plt.axvline(0, c='g')
        plt.axhline(0, c='g')
        for i in range(symbols[self.method].size):
            mask = self.correct == i
            plt.scatter(np.real(self.samples[mask]), np.imag(self.samples[mask]), marker='o', label=str(i), s=25)
        if legend:
            plt.legend()
        plt.xlabel('Real')
        plt.ylabel('Imag')
        plt.title('Constellation')
        plt.grid(True, 'major')
        m = np.max(np.abs(np.hstack((np.real(self.samples), np.imag(self.samples))))) * 1.5
        plt.axis(np.array([-1, 1, -1,  1]) * m)
        plt.show()


class WDMDemodulator(Demodulator):
    def __init__(self, bps=1e6, fs=1e8, filter=[], padding=0):
        super().__init__(bps, fs, filter, padding)

    def setWDM(self, freqCenter=0.0, freqInterval=5e10, channels=1):
        self.freqCenter = freqCenter
        self.freqInterval = freqInterval
        self.channels = channels
        self.freqs = np.linspace(0, self.freqInterval*(self.channels-1), self.channels) - self.freqInterval*(self.channels-1)/2 + self.freqCenter
        return self

    def demodulate(self):
        sup = Demodulator(bps=self.bps, fs=self.fs, filter=self.filter, padding=self.padding)
        self.n = self.correct.shape[-1] + 2*self.padding
        self.Ns = int(self.fs / self.bps)
        self.N = self.Ns * self.n
        self.t = np.arange(0, self.N/self.fs, 1/self.fs)[:self.N]
        self.f = np.arange(-self.fs/2, self.fs/2, self.bps/self.n)[:self.N]

        def ExtractFreq(freqCenter):
            extracted = np.fft.fftshift(np.fft.fft(self.signal * np.exp(-1j * 2 * np.pi * freqCenter * self.t)))
            extracted[np.abs(self.f)>self.freqInterval/2] = 0
            return np.fft.ifft(np.fft.fftshift(extracted))
        
        self.samples = np.zeros(shape=self.correct.shape, dtype=complex)
        self.bits = np.zeros(shape=self.correct.shape, dtype=int)
        for i in range(self.channels):
            sup.setCorrectBits(self.correct[i, :], self.method).setSignal(ExtractFreq(self.freqs[i])).demodulate()
            self.samples[i, :] = sup.samples
            self.bits[i, :] = sup.bits
        return self
    
    def plotConstellation(self, legend = True,
                          codition=lambda x: np.ones(x.shape, dtype=bool)):
        globalMask = codition(self.correct)
        plt.figure(num='Constellation')
        plt.subplot().set_aspect('equal')
        plt.axvline(0, c='g')
        plt.axhline(0, c='g')
        for i in range(symbols[self.method].size):
            mask = (self.correct == i) & globalMask
            plt.scatter(np.real(self.samples[mask]), np.imag(self.samples[mask]), marker='o', label=str(i), s=25)
        if legend:
            plt.legend()
        plt.xlabel('Real')
        plt.ylabel('Imag')
        plt.title('Constellation')
        plt.grid(True, 'major')
        m = np.max(np.abs(np.hstack((np.real(self.samples), np.imag(self.samples))))) * 1.5
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
                  method: Literal['OOK', 'QPSK', '16QAM', '32QAM', '64QAM', 'PAM']='OOK',
                  maxEnergy: float | np.number = None,
                  minEnergy: float | np.number = None,
                  averageEnergy: float | np.number = None,
                  symbolCount: int | np.number = None):
    if method == 'PAM':
        symbols['PAM'] = np.linspace(minEnergy, maxEnergy, symbolCount)
    alphabet = symbols[method]
    if maxEnergy is not None:
        factor = np.sqrt(maxEnergy) / np.max(np.abs(alphabet))
    elif minEnergy is not None:
        if minEnergy != 0:
            factor = np.sqrt(minEnergy) / np.min(np.abs(alphabet))
        else: 
            print('[ERROR] Invalid minEnergy = 0.')
    elif averageEnergy is not None:
        factor = np.sqrt(averageEnergy) / np.mean(np.abs(alphabet))
    else:
        factor = 1.0
    return np.array(alphabet[bits] * factor)