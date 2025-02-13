# Generate waveform based on bits input.

from typing import Literal
import numpy as np
import matplotlib.pyplot as plt
from Demodulator import Demodulator, symbolMapping, WDMDemodulator, PDMDemodulator

# Basic waveform generator
class BasicWaveGenerator:
    # Constructor
    def __init__(self, bps=1e6, fs=1e8, 
                 basicWave:Literal['Square', 'Gaussian', 'RaisedCosine', 'RC']='Gaussian', 
                 alpha:float|Literal['auto']='auto'):
        self.setBps(bps).setFs(fs).setBasicWave(basicWave, alpha).addZero(0)
        self.generated = False

    # Set bps
    def setBps(self, bps):
        self.bps = bps
        self.fs = 20 * bps
        self.generated = False
        return self

    # Set fs
    def setFs(self, fs):
        self.fs = fs
        self.generated = False
        return self
    
    # set bits
    def setBits(self, bits):
        if isinstance(bits, (int, float, np.number)):
            bits = [bits]
        self.raw = np.array(bits)
        self.generated = False
        return self
    
    # set modulation method
    def modulation(self,
                      method: Literal['OOK', 'QPSK', '16QAM', '32QAM', '64QAM', 'PAM']='OOK',
                      maxEnergy: float | np.number=None,
                      minEnergy: float | np.number=None,
                      averageEnergy: float | np.number=None,
                      symbolCount: int | np.number=None):
        self.method = (method, maxEnergy, minEnergy, averageEnergy, symbolCount)
        self.generated = False
        return self
    
    # pad the signal with zeros. Usually 1~4 bits is enough
    def addZero(self, n:int):
        self.padding = n
        self.generated = False
        return self

    def zeroAdder(self, bits, n):
        return np.pad(np.array(bits), (n, n), 'constant', constant_values=0)

    # set waveform type
    def setBasicWave(self, 
                     basicWave:Literal['Square', 'Gaussian', 'RaisedCosine', 'RC']='Gaussian', 
                     alpha:float|Literal['auto']='auto'):
        self.basicWave = basicWave
        if alpha == 'auto':
            if basicWave == 'Gaussian':
                # For Gaussian pulse, alpha = pi * T0, 3dB bandwidth = 0.5887 / alpha
                self.alpha = 0.5887 / self.bps / 2
            elif basicWave == 'RaisedCosine' or basicWave == 'RC':
                self.alpha = 0.5
        else:
            if alpha == 0:
                alpha = 1e-20
            self.alpha = alpha
        return self

    # Generate waveform
    def generate(self, freqType:Literal['all', 'freq', 'omega']='all'):
        self.calcBasicInfo()
        self.bits = symbolMapping(self.raw, *self.method)
        self.bits = self.zeroAdder(self.bits, self.padding)
        self.waveform = self.generator(self.bits)
        self.generated = True
        return self.getWave(freqType)

    # Calc basic information
    def calcBasicInfo(self):
        try:
            self.n = self.raw.shape[-1] + 2*self.padding
            self.Ns = int(self.fs / self.bps)
            self.N = self.Ns * self.n
            self.t = np.arange(0, self.N/self.fs, 1/self.fs)[:self.N]
            self.f = np.arange(-self.fs/2, self.fs/2, self.bps/self.n)[:self.N]
            self.w = 2*np.pi*self.f
        except AttributeError:
            print('Some attributes are missing. Please set bps, fs, and basicWave first.')

    # Basic generator
    def generator(self, bits):
        # Calculate basic impulse
        base = np.zeros(self.Ns)
        base[base.size//2] = 1
        waveform = np.kron(bits, base)
        self.filter = np.ones(self.N)

        # Generate waveform based on basic waveform type
        match self.basicWave:
            case 'Square':
                # Pass the Square filter
                base = np.ones(self.Ns)
                waveform = np.kron(bits, base)
                self.filter = np.fft.fftshift(np.fft.fft(base))
            case 'Gaussian':
                # Pass the Gaussian filter
                self.filter =  np.exp(-self.alpha**2 * self.f**2) * np.sqrt(self.Ns * self.alpha * self.bps * np.sqrt(2 / np.pi))
                waveform = np.fft.ifft(np.fft.fft(waveform) * np.fft.fftshift(self.filter))
            case 'RaisedCosine' | 'RC':
                # Pass the raised cosine filter
                self.filter = np.zeros(self.N)
                upper = np.abs(self.f) <= (1+self.alpha)*self.bps/2
                self.filter[upper] = 0.5*(1 + np.cos(np.pi/self.alpha/self.bps * (np.abs(self.f[upper])-(1-self.alpha)*self.bps/2)))
                self.filter[np.abs(self.f) <= (1-self.alpha)*self.bps] = 1
                waveform = np.fft.ifft(np.fft.fft(waveform) * np.fft.fftshift(self.filter))
            case _:
                print('Unknown basic waveform type')
        return waveform
        
    # Get waveform
    def getWave(self, freqType:Literal['all', 'freq', 'omega']='all'):
        if not self.generated:
            self.generate()
        if freqType == 'all':
            return self.waveform, self.t, self.f, self.w
        elif freqType == 'freq':
            return self.waveform, self.t, self.f
        elif freqType == 'omega':
            return self.waveform, self.t, self.w
        else:
            print('Unknown frequency type')

    # Plot waveform
    def plot(self):
        if not self.generated:
            self.generate()

        plt.figure(num='Waveform')
        plt.subplot(2, 1, 1)
        plt.plot(self.t, np.abs(self.waveform))
        plt.title('Waveform')
        plt.xlabel('Time (s)')
        plt.ylabel('Amplitude')
        plt.grid(True, 'major')

        plt.subplot(2, 1, 2)
        plt.plot(self.f, np.abs(np.fft.fftshift(np.fft.fft(self.waveform))))
        plt.title('Frequency Spectrum')
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('Magnitude')
        plt.grid(True, 'major')
        
        plt.subplots_adjust(hspace=0.5)
        plt.show()
    
    def info(self):
        print(f'bps: {self.bps:.2e}')
        print(f'fs: {self.fs:.2e}Hz')
        print(f'Symbol Period: {1/self.bps*1e12:.2f}ps')
        print('basicWave:', self.basicWave)
        if self.basicWave == 'Gaussian':
            print(f'Half width: {self.alpha / np.pi * 1e12:.2f}ps')
        elif self.basicWave == 'RaisedCosine' or self.basicWave == 'RC':
            print('Alpha: ', self.alpha)
        return self

    def demodulator(self):
        return Demodulator(self.bps, self.fs, self.filter)\
                .setCorrectBits(self.raw, self.method[0])\
                .setPadding(self.padding)


# Wave Generator for WDM
class WDMWaveGenerator(BasicWaveGenerator):
    def __init__(self, bps=1e6, fs=1e8, 
                 basicWave:Literal['Square', 'Gaussian', 'RaisedCosine', 'RC']='Gaussian',
                 alpha:float|Literal['auto']='auto',
                 freqCenter=0.0, freqInterval=5e10, channels=1):
        self.setBps(bps).setFs(fs).setBasicWave(basicWave, alpha).setWDM(freqCenter, freqInterval, channels).addZero(0)
        self.generated = False

    def setWDM(self, freqCenter=0.0, freqInterval=5e10, channels=1):
        self.freqCenter = freqCenter
        self.freqInterval = freqInterval
        self.channels = channels
        self.freqs = np.linspace(0, self.freqInterval*(self.channels-1), self.channels) - self.freqInterval*(self.channels-1)/2 + self.freqCenter
        self.fs = max(self.fs, 2 * np.max(self.freqs) + self.freqInterval) 
        return self
    
    def setBits(self, bits):
        self.raw = np.array(bits)
        self.n = self.raw.shape[1]
        return self
    
    def zeroAdder(self, bits, n):
        return np.pad(np.array(bits), ((0, 0), (n, n)), 'constant', constant_values=0)
    
    def generator(self, bits):
        self.waveform = np.zeros((self.N), dtype=complex)
        for i in range(self.channels):
            self.waveform += super().generator(bits[i, :]) * np.exp(1j * 2*np.pi*self.freqs[i] * self.t)
        return self.waveform
    
    def demodulator(self):
        return WDMDemodulator(self.bps, self.fs, self.filter)\
                .setCorrectBits(self.raw, self.method[0])\
                .setPadding(self.padding)\
                .setWDM(self.freqCenter, self.freqInterval, self.channels)
    
    
# Wave Generator for PDM
class PDMWaveGenerator(WDMWaveGenerator):

    def setBits(self, xbits, ybits=None):
        self.xraw = np.array(xbits)
        if ybits is None:
            self.yraw = np.zeros_like(self.xraw)
        else:
            self.yraw = np.array(ybits)
        return self
    
    def calcBasicInfo(self):
        try:
            self.n = self.xraw.shape[-1] + 2*self.padding
            self.Ns = int(self.fs / self.bps)
            self.N = self.Ns * self.n
            self.t = np.arange(0, self.N/self.fs, 1/self.fs)[:self.N]
            self.f = np.arange(-self.fs/2, self.fs/2, self.bps/self.n)[:self.N]
            self.w = 2*np.pi*self.f
        except AttributeError:
            print('Some attributes are missing. Please set bps, fs, and basicWave first.')

    def generate(self, freqType:Literal['all', 'freq', 'omega']='all'):
        self.calcBasicInfo()
        sup = WDMWaveGenerator(self.bps, self.fs, self.basicWave, self.alpha, self.freqCenter, self.freqInterval, self.channels)
        sup.addZero(self.padding).modulation(*self.method).setBits(self.xraw).generate()
        self.filter = sup.filter
        self.xbits = sup.bits
        xwave = sup.waveform
        sup.setBits(self.yraw).generate()
        self.ybits = sup.bits
        ywave = sup.waveform
        self.waveform = np.vstack((xwave, ywave))
        self.generated = True
        return self.getWave(freqType)
    
    def plot(self):
        if not self.generated:
            self.generate()

        plt.figure(num='X - Waveform')
        plt.subplot(221)
        plt.plot(self.t, np.abs(self.waveform[0, :]))
        plt.title('Waveform')
        plt.xlabel('Time (s)')
        plt.ylabel('Amplitude')
        plt.grid(True, 'major')

        plt.subplot(223)
        plt.plot(self.f, np.abs(np.fft.fftshift(np.fft.fft(self.waveform[0, :]))))
        plt.title('X - Frequency Spectrum')
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('Magnitude')
        plt.grid(True, 'major')

        plt.subplot(222)
        plt.plot(self.t, np.abs(self.waveform[1, :]))
        plt.title('Y - Waveform')
        plt.xlabel('Time (s)')
        plt.ylabel('Amplitude')
        plt.grid(True, 'major')

        plt.subplot(224)
        plt.plot(self.f, np.abs(np.fft.fftshift(np.fft.fft(self.waveform[1, :]))))
        plt.title('Y - Frequency Spectrum')
        plt.xlabel('Frequency (Hz)')
        plt.ylabel('Magnitude')
        plt.grid(True, 'major')
        
        plt.subplots_adjust(hspace=0.5, wspace=0.4)
        plt.show()
    
    def demodulator(self):
        return PDMDemodulator(self.bps, self.fs, self.filter, self.padding)\
                .setWDM(self.freqCenter, self.freqInterval, self.channels)\
                .setCorrectBits(self.xraw, self.yraw, self.method[0])