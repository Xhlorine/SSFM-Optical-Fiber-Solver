from WaveGenerator import BasicWaveGenerator, WDMWaveGenerator
from SSFMSolver import SSFMSolver
from Demodulator import Demodulator

import matplotlib.pyplot as plt
import numpy as np

# Notice:
# 1. When the time span is not wide enough, use function "addZero" to introduce padding bits.
# 2. When the spectrum is not wide enough, raise "fs" to gain higher frequency span.

if __name__ == '__main__':
    # SingleChannel, WDM, PDM
    example = 'WDM'
    if example == 'SingleChannel':
        gen = BasicWaveGenerator(bps=1e9, fs=4e12, basicWave='Gaussian').setBits(np.arange(100))\
                    .addZero(2).modulation('PAM', maxEnergy=100, minEnergy=0, symbolCount=100)
        gen.info()
        signal, t, f, w = gen.generate(freqType='all')
        # gen.plot()
        dem = gen.demodulator()
        dem.setSignal(signal)
        dem.demodulate()
        dem.plotConstellation()

        fiber = SSFMSolver(alphaDB=0, beta=[0, 0], gamma=2, L=100, dz=0.5, title='test').setInput(signal, t, w)
        transmitted = fiber.propagate()
        fiber.plot()

        dem.setSignal(transmitted)
        dem.demodulate()
        dem.plotConstellation()

        plt.figure(num='Phase Rotation')
        plt.subplot(2, 1, 1)
        plt.plot(gen.raw**2, np.unwrap(np.angle(dem.samples)))
        plt.grid(True, 'both')
        plt.subplot(2, 1, 2)
        plt.plot(gen.raw**2, np.abs(dem.samples)**2)
        plt.show()
        print(np.sum(dem.filter**2))
    elif example == 'WDM':
        gen = WDMWaveGenerator(bps=1e9, fs=1e12).setWDM(freqCenter=2.5e10, freqInterval=5e10, channels=3)
        gen.setBits([
            [0, 1, 2, 3],
            [0, 1, 2, 3],
            [0, 1, 2, 3]
        ]).modulation('QPSK', maxEnergy=1e-1).addZero(2).info()
        wave, t, w = gen.generate('omega')

        fiber = SSFMSolver(alphaDB=0, beta=[0, 0], gamma=2, L=80, dz=0.5, title='test').setInput(wave, t, w)
        output = fiber.propagate()
        fiber.plot()

        dem = gen.demodulator()
        dem.setSignal(output)
        dem.demodulate()
        dem.plotConstellation()
    elif example == 'PDM':
        pass
    else:
        print('Invalid example name.')


# TODO:
# 1. How to Demodulate? Refer to the paper from SHJU.
# 2. Test WaveGenerator of WDM and PDM.