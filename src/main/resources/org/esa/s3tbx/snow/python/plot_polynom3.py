import matplotlib.pyplot as plt
import numpy as np

a = 1.938
b = -4.77
c = 7.697
d = -4.092

t = np.arange(0.4, 1.2, 0.01)

s = a + b*t +c*t*t + d*t*t*t
plt.plot(t, s)

data_x = [0.4, 0.753, 0.865, 1.02]
data_y = data_y = [1.0, 0.963, 0.922, 0.737]

plt.plot(data_x, data_y, 'D')

#plt.ylim([0.7, 1.0])
plt.xlabel('Wavelength (microns)')
plt.ylabel('Albedo (dl)')
plt.title('Fit to 3rd order polynom')
plt.grid(True)
plt.savefig("albedo_poly3.png")
plt.show()