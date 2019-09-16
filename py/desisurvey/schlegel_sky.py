import numpy as np
import pylab as pl

from   desisurvey.etc             import schlegel_twi
from   feasiBGS.feasibgs.skymodel import Isky_newKS_twi


if __name__ == '__main__':
  alts    = np.arange(-30., -10., 0.5)
  expfacs = schlegel_twi(alts)

  pl.semilogy(alts, expfacs)

  pl.axhline(2.5, xmin=0., xmax=1., c='k')

  pl.xlim(-30., -10.)
  pl.ylim(0.1,  100.)
  
  pl.show()

  ##
  ##  self._snr2frac = self._snr2frac_start + self.signal ** 2 / self.background / self.texp_total
