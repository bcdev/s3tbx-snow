import matplotlib.pyplot as plt
import numpy as np

a = -2.3
b = 3.3
#lambda_0 = 1326.6
lambda_0 = 1.326
#delta_lambda = 125.3
delta_lambda = 0.1253

#t = np.arange(400.0, 1200.0, 1.0)
t = np.arange(0.4, 1.2, 0.01)

s = a + b / (1.0 + np.exp((t - lambda_0)/delta_lambda))
plt.plot(t, s)

#data_x = [400., 753., 865., 1020.]
data_x = [0.4, 0.753, 0.865, 1.02]
data_y = [1.0, 0.963, 0.922, 0.737]


plt.plot(data_x, data_y, 'D')

plt.ylim([0.7, 1.0])
plt.xlabel('Wavelength (microns)')
plt.ylabel('Albedo (dl)')
plt.title('Fit to AK, eq. (3)')
plt.grid(True)
plt.savefig("albedo_eq3.png")
plt.show()