import numpy as np
import matplotlib.pyplot as plt



# t, y = np.loadtxt("Matlab - per studenti/T20.dat").T


t1, y1 = np.loadtxt("1e01_20C/1atmCO/b1us.dat", skiprows=5, delimiter=",").T
t2, y2 = np.loadtxt("1e01_20C/01atmCO/i1us.dat", skiprows=5, delimiter=",").T

y1 -= y1.mean()
y2 -= y2[:9000].mean()

offset = 6


assert y1[offset:].size == y2[:-offset].size
# plt.plot(y1[:-offset] * 0.8)
plt.plot(y2[offset:])
plt.plot(y2[offset:] - y1[:-offset])
# plt.plot(y2[:-4] - y1[4:])
# plt.plot(y2 + y1)
# plt.plot(t2, y2)
plt.show()