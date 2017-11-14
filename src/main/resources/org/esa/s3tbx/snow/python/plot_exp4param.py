from matplotlib import pyplot
import matplotlib.pyplot as plt
import numpy as np

# f(v; a, kappa1, L, b) := a *(exp(-(2*PI/v) * (kappa1*L + kappa2(v)*L)))^b
a = 0.5968738875747901
kappa_1 = -1.0545344752211868E-6
L = 85750.43560476472
b = 0.3629900776180108

# t = np.arange(0.4, 1.2, 0.01)

# refractive index:
# with open("F:/olafd/s3snow/from_alex/refr_index.csv") as f:
#     datarefr = f.read()
# datarefr = datarefr.split('\n')
#
# im_k = np.asarray([row.split(' ')[1] for row in datarefr]).astype(float)

# s = a + b / (1.0 + np.exp(c*t + d))
# s = a * np.exp(-2.0*np.pi/t * (kappa_1*L + ))
# plt.plot(t, s)

exp4param_x = np.arange(0.4, 1.04, 0.02)
# exp4param_y = [0.94746, 0.9232, 0.90174, 0.88254, 0.86527, 0.8496, 0.8354, 0.8224, 0.8104, 0.7994, 0.7891, 0.7794,
#                0.7704, 0.7618, 0.7536, 0.745, 0.7385, 0.7314, 0.7222, 0.7127,
#                0.7037, 0.6953, 0.6874, 0.6793, 0.6705, 0.6596, 0.6457, 0.6275, 0.6037, 0.5734, 0.5361, 0.4917]
exp4param_y = [.999, .97529,.9537, .93448, .91713, .9014, .88706, .8739, .8618, .8505, .84, .8301,
               .8208, .8119, .8033, .79154, .7882, .7799, .7694, .7584,
               .7479, .7382, .729, .709, .6958, .6784, .6553, .625, 0.6037, 0.5734, 0.5361, 0.4917]
# exp4param_x = [0.4, 0.4125, 0.4425, .49, .51, .56, .62, .665, .67375, .68125, .70875, .75375, .76125,
#                .764375, .7675, .77875, .865, .885, .9, .94, 1.02]
# exp4param_y = [0.94746, 0.9232, 0.89922, 0.85728, 0.84238, 0.8105, 0.7795, 0.75976, 0.75617, 0.75315, 0.742, 0.72523, 0.72168,
#                0.72019, 0.7187, 0.71338, 0.67729, 0.668, 0.65969, 0.6275, 0.4917]

alb_from_brr_vis_x = [0.4, 0.753, 0.865, 1.02]
alb_from_brr_vis_y = [1.0, 0.76812, 0.7226, 0.4834]

# plt.plot(exp4param_x, exp4param_y, 'D')
plt.plot(exp4param_x, exp4param_y)

line2, = plt.plot(alb_from_brr_vis_x, alb_from_brr_vis_y, 'D', label='bla')
data_legend = plt.legend(handles=[line2], loc=1, prop={'size': 10})
data_legend.legendHandles[0].set_linewidth(2.0)
pyplot.gca().add_artist(data_legend)

plt.ylim([0.4, 1.0])
plt.xlabel('Wavelength (microns)')
plt.ylabel('Albedo (dl)')
plt.title('Fit to AK, eq. (3)')
plt.grid(True)
plt.savefig("albedo_eq3.png")
plt.show()