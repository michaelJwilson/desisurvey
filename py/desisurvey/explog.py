import  os 
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
from    desisurvey.desiephem import  Ephemerides
from    desisurvey.utils     import  local_noon_on_date
from    desisurvey.etc       import  bright_exposure_factor
from    datetime             import  datetime, date


year           = 2019
FULLMOONNIGHTS = 7

root           = os.environ['CSCRATCH'] + '/feasiBGS/wilson/'

mayall         = get_location()

##  tfile      = desimodel.io.findfile('footprint/desi-tiles.fits')
tfile          = 'bright-tiles.fits'
tiles          = fits.open(tfile)

tiles          = Table(tiles[1].data)
tiles          = tiles[tiles['PROGRAM'] == 'BRIGHT']
tiles          = tiles[tiles['IN_DESI'] ==        1]

exps           = fits.open(root + 'exposures_surveysim.fits')
exps           = Table(exps[1].data)

exps.sort('MJD')

exps['EXPID']  = np.arange(len(exps['MJD']))

##
exps['KPT']    = (Time([x.value for x in exps['MJD'].quantity], format='mjd') - TimeDelta(60. * 60. * 7.0, format='sec')).iso
exps['DAY']    = np.array([get_date(x).isoformat().split(' ')[0] for x in exps['MJD']])          ##  UTC
exps['YEAR']   = np.array([np.int(x.split('-')[0]) for x in exps['DAY']])                        ##  UTC

##  Year Cut.
isin           = exps['YEAR'] == year
##  exps       = exps[isin]

##  print(exps)

exps         = expderived(exps, tiles=tiles)

exps['TWI']  = np.array([1 if x >= -15. else 0 for x in exps['SUNALT']])

cs           = [SkyCoord(ra = ra * u.degree, dec = dec * u.degree, frame='icrs').transform_to('barycentrictrueecliptic') for ra, dec in zip(exps['RA'].quantity, exps['DEC'].quantity)]
exps['ELON'] = [x.lon.value for x in cs]
exps['ELAT'] = [x.lat.value for x in cs]

exps.write(root + 'exposures_vadd.fits', format='fits', overwrite=True)

##  Low ecliptic latitude cut. 
isin         = [np.abs(x) <= 15. for x in exps['ELAT']]
##exps       = exps[isin]

exps.sort('MJD')

exps.remove_column('SKY')
exps.remove_column('IN_DESI')

## 
dat          = Table(fits.open(root + 'stats_surveysim.fits')[1].data)

dat['DAY']   = np.array([get_date(x).isoformat().split(' ')[0] for x in dat['MJD']])          ##  UTC   
dat['YEAR']  = np.array([np.int(x.split('-')[0]) for x in dat['DAY']])
dat['MID']   = dat['MJD'] + 0.5                                                               ##  MJD at midnight, rather than noon. 

dat['topen'] = np.sum(dat['topen'].quantity, axis=1)
dat['tdead'] = np.sum(dat['tdead'].quantity, axis=1)

##  Populate nights that have no exposures with lunar data.  
##  dat      = expderived(dat, tiles=None, _mjd='MID', unmasked=False)

##  Year cut.
isin         = dat['YEAR'] == year
##  dat      = dat[isin]

dat.remove_column('completed')
dat.remove_column('nexp')
dat.remove_column('nsetup')
dat.remove_column('nsplit')
dat.remove_column('nsetup_abort')
dat.remove_column('nsplit_abort')

##  Load ephem.
##  start        = date(year = 2019, month = 12, day =  1)
##  stop         = date(year = 2024, month = 11, day = 30)
##  fname        = 'ephem_2019-12-01_2024-11-30.fits'

start        = date(year = 2019, month =  1, day =  1)
stop         = date(year = 2025, month = 12, day = 31)

fname        = root + 'ephem_2019-01-01_2025-12-31.fits'

_ephem       = Ephemerides(start, stop, restore=fname)


dat['fmoon'] = np.array([_ephem.is_full_moon(desisurvey.utils.get_date(x), num_nights=FULLMOONNIGHTS) for x in dat['MJD']]).astype(np.int)
dat['msoon'] = np.array([desisurvey.utils.is_monsoon(desisurvey.utils.get_date(x)) for x in dat['MJD']]).astype(np.int)

dat.remove_column('MJD')
dat.remove_column('MID')

##  print(dat)

## 
result      = join(exps, dat, keys='DAY', join_type='inner')

result.remove_column('YEAR_1')
result.remove_column('YEAR_2')

## 
tsciences   = []
tsetups     = []
tsplits     = []

for i in range(len(result['tscience'])):
  try:
    tsciences.append(result['tscience'][i, result['PASS'][i].astype(np.int)])
    tsplits.append(  result['tsplit'][i, result['PASS'][i].astype(np.int)])
    tsetups.append(  result['tsetup'][i, result['PASS'][i].astype(np.int)])
    
  except:
    tsciences.append(np.NaN)
    tsplits.append(np.NaN)
    tsetups.append(np.NaN)
    
result.remove_column('tscience')
result.remove_column('tsplit')
result.remove_column('tsetup')

result['tscience/pass'] = np.array(tsciences)
result['tsplit/pass']   = np.array(tsplits)
result['tsetup/pass']   = np.array(tsetups)

result.sort('MJD')

result.remove_column('MJD')
result.remove_column('STAR_DENSITY')
result.remove_column('PROGRAM')
result.remove_column('OBSCONDITIONS')

##  print(exps)
##  print(dat)
##  print(result)

for i, date in enumerate(dat['DAY']):
  if dat['fmoon'].quantity[i] == 1:
    print('\n\nExposures for day: {}    FULL MOON \n'.format(date))

  else:
    print('\n\nExposures for day: {}\n'.format(date))

  subsample = result[result['DAY'] == date]
  subsample.remove_column('DAY')
  
  print(subsample)
  
  ##  time.sleep(10)

result.write(root + 'results.fits', format='fits', overwrite=True)
    
print('\n\nDone.\n\n')
