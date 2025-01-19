# Recursively Solve NLSE in optical fibers by SSFM
import math
import numpy as np
import matplotlib.pyplot as plt

class SSFMSolver:
    def __init__(self, alphaDB=0.0, beta=[15, 0], gamma=2, L=80, dz=0.5, title=None):
        self.setAlpha(alphaDB).setBeta(beta).setGamma(gamma).setL(L).setDz(dz).usePDM(False).setTitle(title)

    def setAlpha(self, alpha, unit='dB'):
        if unit == 'dB':
            self.alphaDB = alpha
            self.alpha = alpha / 4.343
        elif unit == 'linear' or unit == '1':
            self.alpha = alpha
            self.alphaDB = alpha * 4.343
        else:
            self.alphaDB = alpha
            self.alpha = alpha / 4.343
            print('[WARNING] Undefined unit. Do you mean dB?')
        return self
    
    def setBeta(self, beta=[15, 0]):
        self.beta = np.array(beta)
        return self
    
    def setGamma(self, gamma=2):
        self.gamma = gamma
        return self
    
    def setL(self, L=80):
        self.L = L
        return self
    
    def setDz(self, dz=0.5):
        self.dz = dz
        return self

    def setInput(self, input, t, w):
        self.input = np.asarray(input, dtype=complex)
        self.t = t * 1e12   # convert t from s  to ps
        self.w = w / 1e12   # convert w from Hz to THz
        return self

    def usePDM(self, enable=True):
        self.enablePDM = enable
        return self
    
    def setTitle(self, title):
        self.title = title
        return self
    
    def linearOperator(self, dz):
        dispersion = np.zeros(self.w.shape, dtype=complex)
        for i, b in enumerate(self.beta):
            dispersion +=  1j * (-1)**(i+2) * b / math.factorial(i+2) * self.w**(i+2)
        linear = np.exp((-self.alpha/2 + dispersion) * dz)
        return np.asarray(linear, dtype=complex)

    def nonlinearOperator(self, start):
        return np.asarray(start * np.exp(1j * self.gamma * self.dz * np.abs(start)**2), dtype=complex)
    
    def nonlinearOperatorPDM(self, start):
        return start

    def propagate(self):
        # symmetric SSFM
        if self.enablePDM:
            nonlinearStep = self.nonlinearOperatorPDM
        else:
            nonlinearStep = self.nonlinearOperator
            linearStep = self.linearOperator(self.dz)
            linearStepHalf = self.linearOperator(self.dz / 2)
        self.output = self.input
        # Iterations
        self.output = np.fft.ifft(np.fft.fft(self.output) * np.fft.fftshift(linearStepHalf))
        for i in range(0, int(self.L/self.dz)-1):
            self.output = nonlinearStep(self.output)
            self.output = np.fft.ifft(np.fft.fft(self.output) * np.fft.fftshift(linearStep))
        self.output = nonlinearStep(self.output)
        self.output = np.fft.ifft(np.fft.fft(self.output) * np.fft.fftshift(linearStepHalf))
        # Check parameters
        if np.average(np.abs(self.output[1:10])) > 0.1 * np.max(np.abs(self.output)):
            print('[WARNING] You may need to pad the signal with zeros to avoid edge effects.')
        freq = np.fft.fftshift(np.fft.fft(self.output))
        if np.average(np.max(freq[1:int(freq.size/10)])) > 0.1 * np.max(np.abs(freq)):
            print('[WARNING] You may need to increase fs to avoid aliasing.')
        return self.output
    
    def plot(self):
        if self.title is not None:
            plt.figure(num=self.title)
        else:
            plt.figure()
        plt.subplot(2, 1, 1)
        plt.plot(self.t, np.abs(self.input), label='Input', linestyle='--')
        plt.plot(self.t, np.abs(self.output) * np.exp(self.alpha/2 * self.L), label='Output')
        plt.xlabel('Time (ps)')
        plt.ylabel('Amplitude')
        plt.title('Time Domain')
        plt.grid(True, 'both')
        plt.legend()

        plt.subplot(2, 1, 2)
        plt.plot(self.w/2/np.pi, np.abs(np.fft.fftshift(np.fft.fft(self.input))), label='Input', linestyle='--')
        plt.plot(self.w/2/np.pi, np.abs(np.fft.fftshift(np.fft.fft(self.output * np.exp(self.alpha/2 * self.L)))),
                 label='Output')
        plt.xlabel('Frequency (THz)')
        plt.ylabel('Amplitude')
        plt.title('Frequency Domain')
        plt.grid(True, 'both')
        plt.legend()
        plt.subplots_adjust(hspace=0.5)

        plt.show()

