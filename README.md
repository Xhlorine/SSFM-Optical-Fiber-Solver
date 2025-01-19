# SSFM Optical Fiber Solver

## File List

- `WaveGenerator.py`: Generate waveform based on communication parameters and bits to transmit.
- `SSFMSolver.py`: Simulate the transmission of optical fiber by ***Split Step Fourier Method***.
- `Demodulator.py`: Demodulate and plot constellation.
- `workspace.ipynb`: Detailed introduction in Chinese as well as project "playground".
- `test.py`: A temporary file to test language features for a Python beginner. (Why it's here?)
  > Remember to check your directory before `git add *.py`, please. XD

## Function List

- Generate waveform in the form of single-channel, WDM, PDM.
  > Signal generation of WDM and PDM is not tested yet!
- Transmission of single-channel or WDM signal.
- Coherent demodulation and constellation plot.

## Functions in the Future

- More modulation methods: QPSK, 16QAM, 32QAM, 128QAM.
- More functional interfaces: color map of the whole propagation, animation of constellation rotation, etc.
- Complete mathematical analysis over the Schrödinger Equation and digital simulation.
- Acceleration leveraging NN and other AI tools.


## 省流

光纤传输的数值仿真计算，基于非线性薛定谔方程形式的脉冲传输方程，使用分步傅里叶法计算。

目前可以求解不考虑偏振效应时的光纤传输，并能对单个频带进行解调。

下一步将完善调制、解调方式，并让这几个类用起来更舒适（或许吧，别抱太大希望）。

