import numpy             as     np
import pylab             as     pl
import matplotlib.pyplot as     plt

from   matplotlib.colors import LogNorm
from   plot_brights      import get_brights
from   astropy.table     import Table


print('\n\nWelcome.\n\n')

false_pos, false_neg, true_pos, bright, nbright = get_brights()

cols       = ['PROGRAM', 'MOONFRAC', 'MOONALT', 'SUNALT', 'MIN EXPFAC', 'ISBRIGHT', 'MOONSEP', 'SUNSEP', 'AIRMASS', 'NEXP', 'MJD', 'MOONRA', 'MOONDEC', 'SUNRA', 'SUNDEC', 'LST']

false_pos  = Table(false_pos, names=cols, meta={'name': 'false positive'})
false_neg  = Table(false_neg, names=cols, meta={'name': 'false negative'})
nbright    = Table(nbright,   names=cols, meta={'name': 'not   bright'})

false_pos['LST'][false_pos['LST'] < 0.0] += 360.
false_neg['LST'][false_neg['LST'] < 0.0] += 360.
nbright['LST'][nbright['LST'] < 0.0] += 360.

false_neg.sort('MOONSEP')

##  false_neg  = false_neg[(false_neg['SUNALT'] > -20.) & (false_neg['SUNALT'] < -18.)]

false_nega = false_neg[false_neg['MOONALT'] > 0.]
false_negb = false_neg[false_neg['MOONALT'] < 0.]

print(false_nega)
print('\n\n')
print(false_negb)
print('\n\n')
print(nbright[(nbright['SUNALT'] > -20.) & (nbright['SUNALT'] < -18.)])

##  print(false_pos[false_pos['SUNALT'] < -20.])
##  print('\n\n')
##  print(false_neg[false_neg['SUNALT'] > -20.])

toplot = ['MOONFRAC', 'MOONALT', 'SUNALT', 'MOONSEP', 'SUNSEP']
toplot = ['SUNDEC']
'''
for x in toplot:
    pl.hist(false_nega[x].quantity, bins=20)
    pl.xlabel(x)
    pl.show()

##  bins =  np.arange(-90., 1., 1.0)
'''

diff = false_nega['MOONRA'].quantity - false_nega['LST'].quantity
diff[diff < -180.] += 360.
diff[diff >  180.] -= 360.
diff = np.abs(diff)

plt.scatter(diff, false_nega['MOONRA'].quantity, c=false_nega['MOONSEP'].quantity)
pl.colorbar()
pl.show()

print('\n\nDone.\n\n')
