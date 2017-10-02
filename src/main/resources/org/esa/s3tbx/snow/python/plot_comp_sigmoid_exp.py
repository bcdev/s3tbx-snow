import matplotlib.pyplot as plt
import numpy as np

# 4 data points for following sigmoidal fit:
data_x = [0.4, 0.753, 0.865, 1.02]
data_y = [1.0, 0.963, 0.922, 0.737]
line1, = plt.plot(data_x, data_y, 'D', label='Data at 400, 753, 865, 1020nm')
data_legend = plt.legend(handles=[line1], loc=3, prop={'size': 10})

# sigmoidal approach with 4 fitting parameters:
a = -101521.358793
b = 101522.35987677742
c = 7.58889
d = -20.6012

t = np.arange(0.4, 1.2, 0.01)
s = a + b / (1.0 + np.exp(c * t + d))
plt.plot(t, s)

# exp(sqrt) approach:
t2list = [
    0.400, 0.4125, 0.4425, 0.490, 0.510,
    0.560, 0.620, 0.665, 0.67375, 0.68125,
    0.70875, 0.75375, 0.76125, 0.764375, 0.7675,
    0.77875, 0.865, 0.885, 0.900, 0.940,
    1.02
]
t2 = np.asarray(t2list)
chi_0 = 2.44e-13
lambda_0 = 0.06367
b = 3.62
d = 260.18
s2 = np.exp(-b * np.sqrt(4 * np.pi * d * chi_0 * np.exp(t / lambda_0) / t))

#exp_sqrt_y = [
#0.9962673865827844,0.9959459193168395,0.9950481511844586,0.9931732168906792,0.9921743387911727,
#0.9889583547009275,0.9832387632919509,0.977028319710061,0.9755727977615342,0.974251283164039,
#0.9687586645204144,0.9571226455272059,0.9548005329286031,0.9537964980651493,0.9527703311683053,
#0.9488855637503877,0.906647200550394,0.8928242422137415,0.8812001824950766,0.844153565652906,
#0.7372415128980274
#]
plt.plot(t, s2, 'r')

plt.text(1.0, 0.88, 'Sigmoidal fit:\n'
                    r'$r(\lambda) = \frac{A + B}{1 + \exp(C \lambda + D)}$' '\n'
                    'with' '\n'
                    r'$A = -101521.36$' '\n'
                    r'$B = 101522.36$' '\n'
                    r'$C = 7.59$' '\n'
                    r'$D = -20.6$',
    style='italic', bbox={'facecolor': 'green', 'alpha': 0.5, 'pad': 5})

plt.text(0.5, 0.75, 'Exponential fit:\n'
                    r'$r(\lambda) = \exp(-b \cdot \sqrt {\alpha(\lambda) \cdot d})$' '\n'
                    'with' '\n'
                    r'$\alpha(\lambda) = 4 \pi\chi(\lambda)/\lambda)$' '\n'
                    r'$\chi(\lambda) = \chi_0 exp(\lambda/\lambda_0)$' '\n'
                    r'$\chi_0 = 2.44E-13$' '\n'
                    r'$\lambda_0 = 0.0637 \mu m$' '\n'
                    r'$b = 3.62$' '\n'
                    r'$d = 260.2 \mu m$',
    style='italic', bbox={'facecolor': 'red', 'alpha': 0.5, 'pad': 5})

plt.ylim([0.7, 1.0])
plt.xlabel('Wavelength ($\mu$m)')
plt.ylabel('Albedo (dl)')
plt.title('Comparison of spectral albedo fits')
plt.grid(True)
plt.savefig("albedo_eq3.png")
plt.show()