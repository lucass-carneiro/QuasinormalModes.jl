import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

########################
# Cmd line interaction #
########################

if len(sys.argv) != 2:
    print("Usage: python {} [Number of files to avg.]".format(sys.argv[0]))
    exit(1)
else:
    numFiles = int(sys.argv[1])

###################################
# Load files and compute avg time #
###################################

xData = np.loadtxt("run_1.dat", usecols = 1)
yData = np.zeros(len(xData))

for i in range(1, numFiles + 1, 1):
    yData += 1.0e-9*np.loadtxt("run_%i.dat" % i, usecols = 0)

yData /= numFiles

########
# Plot #
########
font_size = 30
mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['font.family'] = 'Latin Modern Roman'
plt.rcParams['figure.figsize'] = [10, 8]

plt.close('all')
plt.loglog(xData, yData, 'o', color = 'black')

plt.xlabel('Iterations', fontsize = font_size)
plt.ylabel(r'$t$ (s)', fontsize = font_size)

plt.tick_params(axis='both', which='major', labelsize=font_size)

#plt.show()

plt.tight_layout()
plt.savefig('perf.pdf')
