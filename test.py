import numpy as np
import matplotlib.pyplot as plt

# a = np.array([1, 1, 0, 0, 1])
# b = np.ones((1, 6))
# c = np.kron(a, b)
# print(c)

# def GenerateWaveformSquare(bps, fs, bits):
#     # Generate square waveform
#     bits = np.array(bits)
#     n = bits.size
#     t = np.arange(0, bits.size/bps, 1/fs)[:int(bits.size*(fs/bps))]
#     w = np.arange(-1/bps, 1/bps, 2/(fs*bits.size))[:int(bits.size*(fs/bps))]
#     waveform = np.zeros(bits.size)
#     base = np.ones(int(fs/bps))
#     waveform = np.kron(bits, base)
#     waveform = np.sin(2*np.pi*bps*100*t/5)
#     return t, w, waveform

# t, w, wave = GenerateWaveformSquare(1e6, 1e9, [1, 0, 1, 1, 0])

# print(t.size, wave.size)

# plt.figure(num='Waveform')
# plt.plot(t, wave)

# plt.figure()
# plt.plot(w, np.abs(np.fft.fftshift(np.fft.fft(wave))))
# plt.show()

# class Test:
#     def __init__(self, a, b):
#         self.a = a
#         self.b = b

#     def test(self):
#         try:
#             self.c = self.d * 2
#         except AttributeError:
#             print("AttributeError: 'Test' object has no attribute 'd'")

# inst = Test(1, 2)
# inst.test()
# print(inst.a, inst.b, inst.c)

# fs = 1e3
# fc = 100
# N = int(fs)
# t = np.arange(0, 1, 1/fs)[:N]
# f = np.arange(-fs/2, fs/2, fs/N)[:N]

# x = np.sin(2*np.pi*fc*t)

# X = np.fft.fft(x)

# y = np.fft.ifft(X * np.fft.fftshift(np.exp()))

# plt.figure(num='Waveform')
# plt.subplot(2, 1, 1)
# plt.plot(t, x)

# plt.subplot(2, 1, 2)
# plt.plot(t, y)

# plt.show()


# class A:
#     def func(self):
#         print("A.func()")

# class B(A):
#     def func(self):
#         print("B.func()")

# class C(B):
#     def func(self):
#         super(A, self).func()

# inst = C()
# inst.func()


# a = np.array(np.array([[1, 2, 3, 4, 5]]))
# print(a.shape)


# from WaveGenerator import symbolMapping

# syms = np.array([1, 0, 1, 1, 0])
# vals = symbolMapping(syms, 'OOK', maxEnergy=0.1)
# print(vals)


# print(np.linspace(-1, 2, 4))

alphabet = np.array([1, 2, 3, 4])
index = np.array([[0, 1, 2, 3], [3, 2, 1, 0]])
print(index[1, :])