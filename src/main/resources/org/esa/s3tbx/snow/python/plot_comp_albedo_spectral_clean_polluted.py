import matplotlib.pyplot as plt
import numpy as np

# Albedos from file:
with open("F:/olafd/s3snow/from_alex/20181126/comp_albedo_spectral_clean_polluted.csv") as f:
    data = f.read()
data = data.split('\n')

wvl = [row.split(' ')[0] for row in data]
alb_polluted = [row.split(' ')[1] for row in data]
alb_clean = [row.split(' ')[2] for row in data]

# t2 = np.arange(0.4, 2.39, 0.01)
plt.plot(wvl, alb_polluted, 'b')
plt.plot(wvl, alb_clean, 'g')

# plot properties and legends
line1, = plt.plot(wvl, alb_polluted, 'b', label='polluted')
line2, = plt.plot(wvl, alb_clean, 'g', label='clean')
data_legend = plt.legend(handles=[line1,line2], loc=0, prop={'size': 10})
data_legend.legendHandles[0].set_linewidth(2.0)

plt.text(1250.0, 0.6, 'S3A_OL_1_EFR____20180424T100524\n'
                     r''                             '\n'
                     r'                               Lautaret Site' '\n'
                     r'                               45.04N, 6.41E',
         style='italic', bbox={'facecolor': 'green', 'alpha': 0.2, 'pad': 5})

plt.ylim([0.0, 1.0])
plt.xlabel('Wavelength (nm)')
plt.ylabel('Albedo (dl)')
plt.title('Planar spectral albedo from clean and polluted snow algorithm')
plt.grid(True)
plt.savefig("comp_albedo_spectral_clean_polluted.png")

plt.show()
