import pylab              as pl
import numpy              as np
import matplotlib.pyplot  as plt

from   astropy.table      import Table
from   skymodel.sunmodel  import sunmodel
from   skymodel.moonmodel import moonmodel


def get_brights():
    dat       = np.loadtxt('brights.txt')

    ##  by program. 
    darks     = dat[dat[:,0] == 0]
    greys     = dat[dat[:,0] == 1]
    brights   = dat[dat[:,0] == 2]

    ##  bright by etc. 
    bright    = dat[dat[:,5] == 1]
    nbright   = dat[dat[:,5] == 0]

    ##  not bright, but program 2. 
    false_pos = nbright[nbright[:,0] == 2]  ##  Pass 2, not bright. 
    false_neg =  bright[ bright[:,0] != 2]  ##  Not Pass 2, bright. 
    
    true_pos  =  bright[ bright[:,0] == 2]  ##  Pass 2, bright.  
    true_neg  =  bright[ bright[:,0] == 2]  ##  Not Pass2, not bright. i.e. other.   

    return  false_pos, false_neg, true_pos, bright, nbright

def print_tables(dat, names):
    dat       = Table(dat, names=names)
    dat.sort('MIN EXPFAC')

    print
    print
    print(dat)
    
    
if __name__ == '__main__':
    names     = ['PROGRAM', 'MOONFRAC', 'MOONALT', 'SUNALT', 'MIN EXPFAC', 'ISBRIGHT', 'MOONSEP', 'SUNSEP', 'AIRMASS', 'EFFTIME', 'MJD', 'MOONRA', 'MOONDEC', 'SUNRA', 'SUNDEC', 'LST'] 

    false_pos, false_neg, true_pos, bright, nbright = get_brights()

    print_tables(false_pos, names)
    print_tables(false_neg, names)
    
    ##
    nexp_fp   =  np.sum(false_pos[:,9])
    nexp_fn   =  np.sum(false_neg[:,9])
    nexp_tp   =  np.sum( true_pos[:,9])
    
    print('\n\nFalse +, False -, True +.')
    print(nexp_fp, nexp_fn, nexp_tp)
    
    ##
    plt.figure(figsize=(15,7), dpi=200)
        
    pl.plot(false_pos[:,2], false_pos[:,1], '.', markersize=3, label='False + ({0:.1f})'.format(nexp_fp), alpha=0.7)
    pl.plot(false_neg[:,2], false_neg[:,1], '.', markersize=3, label='False - ({0:.1f})'.format(nexp_fn), alpha=0.7)
    pl.plot(true_pos[:,2],  true_pos[:,1],  '.', markersize=3, label='True  + ({0:.1f})'.format(nexp_tp), alpha=0.7)

    ##
    abscissae  = np.arange(0., 90., 0.1)
    ordinates  = 0.54 + 3. / (abscissae + 10.)
    pl.plot(abscissae[ordinates > 0.6], ordinates[ordinates > 0.6], 'k--')

    ## 
    _ordinates = 25. / abscissae
    pl.plot(abscissae[_ordinates <= ordinates], _ordinates[_ordinates <= ordinates], 'k--')

    _ordinates = 30. / abscissae
    pl.plot(abscissae[_ordinates <= 0.6], _ordinates[_ordinates <= 0.6], 'k', alpha=0.6)

    pl.axhline(0.6, xmin=0, xmax=.555, c='k', alpha=0.6)

    pl.xlim(0., 90.)
    pl.ylim(0., 1.0)

    pl.xlabel('Moon altitude.')
    pl.ylabel('Moon illumination.')

    pl.legend(loc=4, frameon=False)

    pl.show()
    ##  pl.savefig('brights.pdf', dpi=800)
