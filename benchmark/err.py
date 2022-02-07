import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

##################
# Load data file #
##################
data_re = np.loadtxt("run_1.dat", usecols = (1,9))

xData_re = data_re[:,0]
yData_re = abs(data_re[:,1])

data_im = np.loadtxt("run_1.dat", usecols = (1,10))

xData_im = data_im[:,0]
yData_im = abs(data_im[:,1])

########
# Plot #
########
font_size = 30
mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['font.family'] = 'Latin Modern Roman'
plt.rcParams['figure.figsize'] = [10, 8]

plt.close('all')
plt.semilogy(xData_re, yData_re, 'o', color = 'black', label = r'$\Re(\omega)$')
plt.semilogy(xData_im, yData_im, 'x', color = 'red', label = r'$\Im(\omega)$')

plt.legend(loc = 'upper right', fontsize = font_size)

plt.xlabel('Iterations', fontsize = font_size)
plt.ylabel(r'$\log|\varepsilon|$', fontsize = font_size)

plt.tick_params(axis='both', which='major', labelsize=font_size)

#plt.show()

plt.tight_layout()
plt.savefig('err.pdf')

