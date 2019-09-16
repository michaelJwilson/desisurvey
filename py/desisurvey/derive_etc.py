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
from    datetime             import  datetime, date


compute = True

root    = os.environ['CSCRATCH'] + '/feasiBGS/wilson/'

tiles   = desisurvey.tiles.get_tiles('./bright-tiles.fits')
result  = Table(fits.open(root + 'results.fits')[1].data)

if compute:
  result['BTIME']       = np.zeros_like(result['TILEID'].quantity, dtype=np.float32)

  for i, tile in enumerate(result['TILEID']):
    result['BTIME'][i]  = 150.  ##  Nominal dark-time exposure.

    result['BTIME'][i] *= desisurvey.etc.bright_exposure_factor(np.array(result['MOONFRAC'].quantity[i].value), np.array(result['MOONALT'].quantity[i].value),\
                                                                np.array(result['MOONSEP'].quantity[i].value), result['SUNALT'].quantity[i].value,\
                                                                np.array(result['SUNSEP'].quantity[i].value), np.array(result['AIRMASS'].quantity[i].value))
    
  ##  Rederive desired exposure time.
  result['ETCTIME']      = result['BTIME']  ##  Bright factor and nominal time.
  result['ETCTIME']     *= desisurvey.etc.dust_exposure_factor(result['EBV_MED'].quantity.value)
  ##  result['ETCTIME'] *= desisurvey.etc.airmass_exposure_factor(result['AIRMASS'].quantity.value)
  ##  result['ETCTIME'] /= ETC.weather_factor(result['SEEING'].quantity.value, result['TRANSP'].quantity.value)

  result['OPTIME'  ]     = np.zeros_like(result['ETCTIME'].quantity)
  result['NEXP']         = np.zeros_like(result['ETCTIME'].quantity)

  exptiles, _            = np.unique(result['TILEID'].quantity, return_counts=True)
  
  for i, tile in enumerate(exptiles):
    isin                   = result['TILEID'].quantity == tile
    
    result['OPTIME'][isin] = np.sum(result['EXPTIME'][isin])
    result['NEXP'][isin]   = len(result[isin])

  result.write(root + 'results_vadd_etc.fits', format='fits', overwrite=True)
  
else:
  result  = Table(fits.open(root + 'results_vadd_etc.fits')[1].data)

##
result.sort('EXPID')

result['EXPTIME'][:] /= 60.
result['ETCTIME'][:] /= 60.
result['BTIME'][:]   /= 60.
result['OPTIME'][:]  /= 60.

result.remove_column('DAY')
result.remove_column('ELON')
result.remove_column('EXPID')
result.remove_column('SUNRA')
result.remove_column('SUNDEC')
result.remove_column('topen')
result.remove_column('tdead')
result.remove_column('tsched')
result.remove_column('tscience/pass')
result.remove_column('tsplit/pass')
result.remove_column('tsetup/pass')
result.remove_column('EXPOSEFAC')
result.remove_column('EBV_MED')
result.remove_column('ELAT')
result.remove_column('IN_IMAGING')
result.remove_column('IMAGEFRAC_GRZ')
result.remove_column('IMAGEFRAC_GR')
result.remove_column('IMAGEFRAC_G')
result.remove_column('IMAGEFRAC_R')
result.remove_column('IMAGEFRAC_Z')
result.remove_column('fmoon')
result.remove_column('msoon')
result.remove_column('BRIGHTVTMAG')
result.remove_column('BRIGHTDEC')
result.remove_column('BRIGHTRA')
result.remove_column('CENTERID')
result.remove_column('RA')
result.remove_column('DEC')

result.sort('KPT')

print(result[result['TILEID'] ==  32193])

result.sort('OPTIME')

##
# result = result[np.abs(result['ELAT'].quantity) < 10.0]
result = result[result['SNR2FRAC'] > 1.]

print(result)

pl.hist(result['OPTIME'], bins=20)
pl.show()

exit(1)

##  Total time:  Total time the shutter was open on this exposure.
passes       = np.unique(result['PASS'].quantity,   return_counts=True)
exptiles, _  = np.unique(result['TILEID'].quantity, return_counts=True)

etc          = np.zeros_like(exptiles, dtype=np.float32)
snr          = np.zeros_like(exptiles, dtype=np.float32)
totexp       = np.zeros_like(exptiles, dtype=np.float32)

for i, tile in enumerate(exptiles):
  isin                   = result['TILEID'].quantity == tile
  
  totexp[i]              = np.sum(result['EXPTIME'][isin])
  snr[i]                 = np.max(result['SNR2FRAC'][isin])     
  '''
  etc[i]                 = 150.
  etc[i]                *= 1.33
  etc[i]                *= desisurvey.etc.dust_exposure_factor(result['EBV_MED'].quantity[result['TILEID']  == tile][0])                                                                                                                                                                                               
  etc[i]                *= desisurvey.etc.airmass_exposure_factor(result['AIRMASS'].quantity[result['TILEID'] == tile][0])
  etc[i]                /= 60.  ##  mins.
  '''
  etc[i]                 = np.float32(result['ETCTIME'][isin][0])                                                                                                                                                                                                                                                                                                                                                                         
  #if np.abs(totexp[i] / etc[i]) < .2:                                                                                                                                                                                                                                                                                                                                                                                               
  print()
  print()
  print(result[isin])

  result['ETCTIME'][isin] = etc[i]

##
plt.figure(figsize=(20, 10))

pl.plot(np.arange(0., 1.e3, 1.), np.arange(0., 1.e3, 1.), alpha=0.3, c='k')
plt.scatter(etc, totexp, c=snr, s=1) ##  vmin=0.6, vmax=1.0

pl.xlim(3., 3.e1)
pl.ylim(3., 3.e1)

pl.colorbar()

pl.xlabel('ETC time [min.]')
pl.ylabel('Open shutter time [min.]')

pl.show()
