import matplotlib.pyplot as plt
import numpy as np

# sigmoidal fit with 4 fitting parameters:
a = -101521.358793
b = 101522.35987677742
c = 7.58889
d = -20.6012

t = np.arange(0.4, 1.2, 0.01)
s = a + b / (1.0 + np.exp(c * t + d))
plt.plot(t, s, 'g')

# exponential expression erivedfrom grain diameter:
chi_0 = 2.44e-13
lambda_0 = 0.06367
b = 3.62
d = 260.18
s2 = np.exp(-b * np.sqrt(4 * np.pi * d * chi_0 * np.exp(t / lambda_0) / t))
plt.plot(t, s2, 'r')

plt.text(0.75, 0.6, 'Sigmoidal fit:\n'
                    r'$r(\lambda) = \frac{A + B}{1 + \exp(C \lambda + D)}$' '\n'
                    'with' '\n'
                    r'$A = -101521.36$' '\n'
                    r'$B = 101522.36$' '\n'
                    r'$C = 7.59$' '\n'
                    r'$D = -20.6$',
    style='italic', bbox={'facecolor': 'green', 'alpha': 0.5, 'pad': 5})

plt.text(0.45, 0.6, 'Exponential fit:\n'
                    r'$r(\lambda) = \exp(-b \cdot \sqrt {\alpha(\lambda) \cdot d})$' '\n'
                    'with' '\n'
                    r'$\alpha(\lambda) = 4 \pi\chi(\lambda)/\lambda)$' '\n'
                    r'$\chi(\lambda) = \chi_0 exp(\lambda/\lambda_0)$' '\n'
                    r'$\chi_0 = 2.44E-13$' '\n'
                    r'$\lambda_0 = 0.0637 \mu m$' '\n'
                    r'$b = 3.62$' '\n'
                    r'$d = 260.2 \mu m$',
    style='italic', bbox={'facecolor': 'red', 'alpha': 0.5, 'pad': 5})

plt.ylim([0.5, 1.0])
plt.xlabel('Wavelength ($\mu$m)')
plt.ylabel('Albedo (dl)')
plt.title('Comparison of spectral albedo fits')
plt.grid(True)
plt.savefig("comp_sigmoid_exp.png")

plt.show()
