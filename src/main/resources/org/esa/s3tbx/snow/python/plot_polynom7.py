import matplotlib.pyplot as plt
import numpy as np

a = 1.03
b = 0.
c = -0.2515
d = 0.
e = 0.4111
f = 0.
g = 0.
h = -0.41505

t = np.arange(0.4, 1.2, 0.01)

s = a + b*t +c*t*t + d*t*t*t + e*t*t*t*t + f*t*t*t*t*t + g*t*t*t*t*t*t + h*t*t*t*t*t*t*t
plt.plot(t, s)

data_x = [0.4, 0.753, 0.865, 1.02]
data_y = [1.0, 0.963, 0.922, 0.737]

plt.plot(data_x, data_y, 'D')

plt.ylim([0.7, 1.0])
plt.xlabel('Wavelength (microns)')
plt.ylabel('Albedo (dl)')
plt.title('Fit to 7th order polynom')
plt.grid(True)
plt.savefig("albedo_poly7.png")
plt.show()