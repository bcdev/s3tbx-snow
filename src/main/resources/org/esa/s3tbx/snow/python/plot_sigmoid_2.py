import matplotlib.pyplot as plt
import numpy as np

a = -101521.358793
b = 101522.35987677742
c = 7.58889
d = -20.6012

# approach with 2 fitting parameters seems sufficient:
a = 0.0
b = 1.0
c = 8.94416
d = -10.15649

t = np.arange(0.4, 1.2, 0.01)

s = a + b / (1.0 + np.exp(c*t + d))
plt.plot(t, s)

data_x = [0.4, 0.753, 0.865, 1.02]
data_y = [1.0, 0.963, 0.922, 0.737]

plt.plot(data_x, data_y, 'D')

#plt.text(0.5, 0.85, 'A = -101521.36\nB = 101522.36\nC = 7.59\nD = -20.60', style='italic',
#    bbox={'facecolor':'red', 'alpha':0.5, 'pad':5})
plt.text(0.5, 0.85, 'A = 8.95\nB = -10.15', style='italic',
    bbox={'facecolor':'red', 'alpha':0.5, 'pad':5})

plt.ylim([0.7, 1.0])
plt.xlabel('Wavelength (microns)')
plt.ylabel('Albedo (dl)')
#plt.title('Fit to Sigmoidal: y = A + B / (1 + exp(Cx+D))')
plt.title('Fit to Sigmoidal: y = 1 / (1 + exp(Ax+B))')
plt.grid(True)
plt.savefig("albedo_eq3.png")
plt.show()