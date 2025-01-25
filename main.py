from WaveGenerator import BasicWaveGenerator
from SSFMSolver import SSFMSolver
from Demodulator import Demodulator

import matplotlib.pyplot as plt
import numpy as np

# Notice:
# 1. When the time span is not wide enough, use function "addZero" to introduce padding bits.
# 2. When the spectrum is not wide enough, raise "fs" to gain higher frequency span.

if __name__ == '__main__':
    # SingleChannel, WDM, PDM
    example = 'SingleChannel'
    if example == 'SingleChannel':
        gen = BasicWaveGenerator(bps=1e9, fs=4e12, basicWave='Gaussian').setBits([0, 0, 1, 1, 2, 2, 3, 3, 2, 2, 1, 2, 0])\
                    .addZero(2).modulation('QPSK', maxEnergy=2)
        gen.info()
        signal, t, f, w = gen.generate(freqType='all')
        # gen.plot()
        dem = gen.demodulator()
        dem.setSignal(signal)
        dem.demodulate()
        dem.plotConstellation()

        fiber = SSFMSolver(alphaDB=0.2, beta=[3, 0], gamma=2, L=100, dz=0.5, title='test').setInput(signal, t, w)
        fiber.propagate()

        dem.setSignal(fiber.output)
        dem.demodulate()
        dem.plotConstellation()
    elif example == 'WDM':
        pass
    elif example == 'PDM':
        pass
    else:
        print('Invalid example name.')


# TODO:
# 1. How to Demodulate? Refer to the paper from SHJU.
# 2. Test WaveGenerator of WDM and PDM.
# 3. Add property "padding" to de-addZero.