from matplotlib import pyplot
import matplotlib.pyplot as plt
import numpy as np

# 4th order spline fit, 400-525nm:
a = 4.913273573521765E-8
b = -5.446912137768008E-7
c = 2.2340617964072435E-6
d = -4.029955265312061E-6
e = 2.70443404633965E-6

t = np.arange(0.4, 0.525, 0.001)
s = a + b*t +c*t*t + d*t*t*t + e*t*t*t*t
plt.semilogy(t, s, 'r', linestyle='-', linewidth=4.0)

# 4th order spline fit, 525-700nm:
a = 7.871975805631149E-7
b = -5.5075218680064285E-6
c = 1.4735203245000926E-5
d = -1.7996480996799473E-5
e = 8.533638647997086E-6

t = np.arange(0.525, 0.71, 0.001)
s = a + b*t +c*t*t + d*t*t*t + e*t*t*t*t
plt.semilogy(t, s, 'r', linestyle='-', linewidth=4.0)

# 4th order spline fit, 700-1020nm:
a = 3.20662332955914E-4
b = -0.0016037321152243488
c = 0.0030016823137293783
d = -0.0024933767575023914
e = 7.763358059437843E-4

t = np.arange(0.71, 1.02, 0.001)
s = a + b*t +c*t*t + d*t*t*t + e*t*t*t*t
line1, = plt.semilogy(t, s, 'r', linestyle='-', linewidth=4.0, label='Piecewise 4th order spline fit')
data_legend = plt.legend(handles=[line1], loc=2, bbox_to_anchor=(0.03, 0.92), prop={'size': 10})
data_legend.legendHandles[0].set_linewidth(2.0)
pyplot.gca().add_artist(data_legend)

# https://atmos.washington.edu/ice_optical_constants data:
with open("F:/olafd/s3snow/from_alex/refr_index.csv") as f:
    data = f.read()
data = data.split('\n')

im_k = np.asarray([row.split(' ')[1] for row in data]).astype(float)

t = np.arange(0.4, 1.04, 0.01)
line2, = plt.semilogy(t, im_k, 'D', label='https://atmos.washington.edu/ice_optical_constants/IOP_2008_ASCIItable.dat')
data_legend = plt.legend(handles=[line2], loc=1, prop={'size': 10})
data_legend.legendHandles[0].set_linewidth(2.0)
pyplot.gca().add_artist(data_legend)

plt.ylim([1.E-11, 1.E-5])
plt.xlabel('Wavelength ($\mu$m)')
plt.ylabel('Refractive Index $\kappa_{2}(\lambda)$')
plt.title('Refractive index fitting')
plt.grid(True)
plt.savefig("refr_index_spline_fitting.png")

plt.show()
