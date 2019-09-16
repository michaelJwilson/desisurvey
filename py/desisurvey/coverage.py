import  time
import  desisurvey
import  astropy.io.fits      as      fits
import  surveysim.stats      as      simstats
import  numpy                as      np
import  astropy.units        as      u
import  pylab                as      pl
import  matplotlib.pyplot    as      plt
import  desimodel

from    astropy.table        import  Table, join
from    astropy.time         import  Time, TimeDelta
from    datetime             import  datetime
from    desisurvey.utils     import  get_location, get_airmass, local_noon_on_date, get_date
from    expderived           import  expderived
from    astropy.coordinates  import  SkyCoord, BaseEclipticFrame
from    desisurvey.ephem     import  Ephemerides
from    desisurvey.utils     import  local_noon_on_date
from    datetime             import  datetime, date


tfile             = desimodel.io.findfile('footprint/desi-tiles.fits')
tiles             = fits.open(tfile)

tiles             = Table(tiles[1].data)
tiles             = tiles[tiles['PROGRAM'] == 'BRIGHT']

##  NOTE:  not all tiles are in DESI in greater passes. 
tiles['CENTERID'] = tiles['TILEID'] - tiles['TILEID'].min()
tiles['CENTERID'] = tiles['CENTERID'] %	5762  ##  TILES PER PASS. 

## 
dat              = fits.open('results_vadd_etc.fits')
dat              = Table(dat[1].data)

dat.remove_column('tsched')
dat.remove_column('topen')
dat.remove_column('tdead')
dat.remove_column('ELON')
dat.remove_column('DAY')
dat.remove_column('tsplit/pass')
dat.remove_column('tscience/pass')
dat.remove_column('tsetup/pass')
dat.remove_column('EXPOSEFAC')

##  NOTE:  not all tiles are in DESI in greater passes.
dat['CENTERID']  = dat['TILEID'] - tiles['TILEID'].min()
dat['CENTERID']  = dat['CENTERID'] % 5762  ##  TILES PER PASS.     

passes           = np.unique(dat['PASS'].quantity,   return_counts=True)
exptiles, _      = np.unique(dat['TILEID'].quantity, return_counts=True)

iscomplete       = dat[dat['SNR2FRAC'].quantity >= 1.0]
cids, cnts       = np.unique(iscomplete['CENTERID'].quantity, return_counts=True)

print('Percentage of three pass coverage: {}'.format(100. * np.sum(cnts == 3) / 2042))
print('Percentage of two pass coverage: {}'.format(100. * np.sum(cnts == 2) / 2042))
print('Percentage of one pass coverage: {}'.format(100. * np.sum(cnts == 1) / 2042))
print('Percentage of zero three pass coverage: {}\n\n'.format(100. * np.sum(cnts == 0) / 2042))

plt.figure(figsize=(20, 10))

##                                                                                                                                                                                                                                                                                                                        
isin = (dat['TILEID'].quantity == dat['TILEID'][np.argmax(dat['OPTIME'].quantity)])

print(dat[isin])
print('\n\n')

##  pl.plot(dat[isin]['RA'] - 90., dat[isin]['DEC'], 'm*', markersize=5)

##  Limit to lowest pass. 
dat              = dat[dat['PASS'].quantity == 5]

##  print(dat)

pl.plot(dat['MOONRA'].quantity, dat['MOONDEC'].quantity, 'ko', markersize=1, alpha=0.4)

np.savetxt('ecliptic-curve.txt', np.c_[dat['MOONRA'].quantity, dat['MOONDEC'].quantity], fmt='%.4lf')

##  plt.scatter(dat['RA'].quantity, dat['DEC'].quantity, c=np.log10(dat['ETCTIME']), marker='o', s=5, cmap='magma')
##  plt.scatter(dat['RA'].quantity, dat['DEC'].quantity, c = (dat['ETCTIME'] / dat['OPTIME']), marker='o', s=5, cmap='magma')    

for coverage, color in zip([1, 2, 3], ['r', 'y', 'g']):
  isin           = [x in cids[cnts == coverage] for x in dat['CENTERID'].quantity]

  ##  print(dat[isin])

  sample         = dat[isin]

  print('\n\n')
  print(sample)

  sample['RA'] -= 90.
  sample['RA'][sample['RA'] < 0.] += 360.
  
  pl.plot(sample['RA'].quantity, sample['DEC'].quantity, color + 'o', markersize=5)

##  pl.colorbar()

print('\n\nDone.\n\n')

pl.show()
