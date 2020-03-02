import  copy
import  ephem
import  desisurvey
import  matplotlib                 as      mpl 
import  astropy.io.fits            as      fits
import  pylab                      as      pl
import  numpy                      as      np
import  astropy.units              as      u
import  matplotlib.pyplot          as      plt
import  matplotlib.cm              as      cm

from    astropy.table              import  Table

from    desisurvey.desiephem       import  Ephemerides
from    desisurvey.utils           import  local_noon_on_date, get_date
from    datetime                   import  datetime, date
from    astropy.time               import  Time
from    astropy.coordinates        import  SkyCoord, EarthLocation, AltAz
from    desitarget.geomask         import  circles
from    desisurvey                 import  plan
from    desisurvey.moonmodel_7160  import  moonmodel
from    desisurvey.sunmodel_7160   import  sunmodel


def whatprogram(mjd, programs, changes):
  index = 0 

  while mjd >= changes[index + 1]:
    index += 1
     
  return  programs[index]

##  prefix       = 'op'
prefix           = 'y1'
##  prefix       = 'y5'

verbose          = True

##  Dark time in nominal conditions.
fidexp           = 150.

##  Get Eddie's tiles -> remppaed to a format similar to the old tiles, e.g. in pass ordering.
##  Center ID defined for pass 5 (Gray) in this instance.  
##  tiles        = Table(fits.open('/global/cscratch1/sd/mjwilson/svdc2019c2/survey/tiles/schlafly/opall.fits')[1].data)

tiles            = Table(fits.open('bright-tiles.fits')[1].data)
tiles            = tiles[tiles['IN_DESI'].quantity > 0]
tiles            = tiles[tiles['PASS'].quantity ==  0]

tiles.sort('CENTERID')

##  print(tiles)

##
cids             = np.unique(tiles['CENTERID'].quantity)

##  Write.                                                                                                                                                                               
np.savetxt('visibility-{}/cids.txt'.format(prefix), cids, fmt='%d')

##  design_hrangle = plan.load_design_hourangle()               

##  config derived constraints.                                                                                                                                                         
config           = desisurvey.config.Configuration(file_name='/Users/MJWilson/Work/desi/repos/chh/desisurvey/py/desisurvey/mjw-config.yaml')
full_moon_nights = config.full_moon_nights()

##  Config. defined.
##  first_day    = config.first_day().isoformat()
##  last_day     = config.last_day().isoformat()

##  One-year rewrite.                                                                                                                                                                      
first_day        = '2020-01-01'
last_day         = '2021-01-01'  

##  Five-year rewrite.
##  first_day    = '2020-06-01'   
##  last_day     = '2025-12-31'

print(first_day, last_day)

first_mjd        = Time(first_day, format='iso').mjd
last_mjd         = Time(last_day,  format='iso').mjd

##  
min_altitude     = config.min_altitude().value

lat              = config.location.latitude()
lon              = config.location.longitude()
elv              = config.location.elevation()

print(lat, lon, elv)

avoid_bodies     = {}

##  Preserve order. 
bodies           = list(config.avoid_bodies.keys)

for body in bodies:
  avoid_bodies[body] = getattr(config.avoid_bodies, body)().to(u.deg)

##  Planet exclusion evaluated at midnight, moon exclusion at each mjd. 
bodies.remove('moon')
  
##
##  mayall        = EarthLocation(lat=lat, lon=lon, height=elv)
mayall            = EarthLocation.from_geodetic(lat=lat, lon=lon, height=elv)

emayall           = ephem.Observer()
emayall.lat       = lat.to(u.rad).value
emayall.lon       = lon.to(u.rad).value
emayall.elevation = elv.to(u.m).value

##
ra                = tiles['RA'].quantity
dec               = tiles['DEC'].quantity

# Calculate the maximum |HA| in degrees allowed for each tile to stay above the survey minimum altitude (plus a 1 deg padding).
cosz_min          = np.cos(90. * u.deg - (config.min_altitude() + 1. * u.deg))
cosHA_min         = ((cosz_min - np.sin(dec * u.deg) * np.sin(lat)) / (np.cos(dec * u.deg) * np.cos(lat))).value
max_abs_ha        = np.degrees(np.arccos(cosHA_min))

print(min_altitude, max_abs_ha.min() / 15., max_abs_ha.max() / 15.)

##  ephem table duration.
start            = date(year = 2019, month =  1,  day = 1)
stop             = date(year = 2025, month = 12, day = 31)

##  Load ephem. 
dat              = Ephemerides(start, stop, restore='/global/cscratch1/sd/mjwilson/feasiBGS/wilson/ephem_2019-01-01_2025-12-31.fits')

##  print(dat._table.columns)
##  print(dat._table)

##  Survey sim derived program hours for each night. 
hrs              = dat.get_program_hours(include_twilight=False)

print('\n\nProgram hour breakdown.')

for i, _hrs in enumerate(hrs.T):
  print('Night  {0:4}:    {1:.4f} Dark;    {2:.4f}  Gray;    {3:.4f}  Bright.'.format(i, _hrs[0], _hrs[1], _hrs[2]))
 
##  Choose same times as those solved for in ephem, but more finely sample than 1/hr due to twilight changes.   
N                = 96 * 3
dt               = 24. / N
t_obj            = np.linspace(0., 1., N + 1)

print('\n\nHours visible for each CENTERID, program and noon-to-noon day.  {} min. sampling.\n\n'.format(dt * 60.))

##  Initialise tracking.
nnights          = 0

##  For each tile, and each night, record the hrs visible in each program. 
hrs_visible      = np.zeros(3 * len(tiles['RA'].quantity) * len(dat._table['noon'].quantity), dtype=np.float).reshape(len(dat._table['noon'].quantity), len(tiles['RA'].quantity), 3)
moons            = []

isonoons         = []
brights          = []

output           = []

for i, noon in enumerate(dat._table['noon'].quantity):
 isonoon     = Time(noon.value, format='mjd').iso.split(' ')[0]

 fullmoon    = dat.is_full_moon(get_date(isonoon))         
 monsoon     = desisurvey.utils.is_monsoon(get_date(isonoon))  

 if (isonoon >= first_day) & (isonoon <= last_day):
  if fullmoon:
    print('{} FULLMOON'.format(isonoon))
    continue

  if monsoon:
    print('{} MONSOON'.format(isonoon))
    continue

  ##  Count the nights that made it. 
  nnights   += 1
  
  ##  This ephem row made the cut. 
  output.append(i)

  ##  Isonoon for this row. 
  isonoons.append([i, isonoon, nnights])

  ##  MJD for midnight on this date. 
  midnight     = noon.value + 0.5
  
  dusk         = dat._table['dusk'].quantity[i]
  dawn         = dat._table['dawn'].quantity[i]

  ##  Includes twilight.                                                                                                                                                                                                              
  bdusk        = dat._table['brightdusk'].quantity[i]
  bdawn        = dat._table['brightdawn'].quantity[i]
  
  ##  Assume this does not vary much during the night.
  ##  moonfrac = dat._table['moon_illum_frac'].quantity[i]  
  
  ##  Calculate LST and hour angle.
  MJD0, MJD1   = bdusk, bdawn
  LST0, LST1   = dat._table['brightdusk_LST'][i], dat._table['brightdawn_LST'][i]
  dLST         = (LST1 - LST0) / (MJD1 - MJD0)

  ##  programs, changes = dat._table['programs'].quantity[i],  dat._table['changes'].quantity[i]                                                                                                                                      
  programs, changes = dat.get_night_program(get_date(isonoon), include_twilight=False, program_as_int=True)  ##  Call bookended by (b)dusk, (b)dawn at either end.    
  
  ##  Sorted by definition.
  indices           = list(range(len(ra)))
  
  ##  Planet exclusion evaluated at midnight, moon exclusion at each mjd. 
  ##  Note:  no solar exclusion currently.
  for body in bodies:
    bdec, bra  = desisurvey.desiephem.get_object_interpolator(dat._table[i], body, altaz=False)(midnight)
    too_close  = desisurvey.utils.separation_matrix([bra] * u.deg, [bdec] * u.deg, ra[indices] * u.deg, dec[indices] * u.deg, avoid_bodies[body])[0]
    
    toremove   = np.argwhere(too_close == True)[:,0]

    ##  Start at the back, otherwise original indexing fails after elements have been removed. 
    for index in sorted(toremove, reverse=True):
      indices.pop(index)
  
  ##  On this day, enumerate over the time samples and determine which CENTERIDs are visible at each time. 
  for j, t in enumerate(t_obj):
    mjd          = noon.value + t

    ##  Include twilight with bdusk and bdawn. 
    if (mjd < dusk) or (mjd > dawn):      
      continue

    program      = whatprogram(mjd, programs, changes)    

    time         = Time(mjd, format='mjd')  ##  UTC.

    moonfrac     = np.float(dat.get_moon_illuminated_fraction(mjd))
    
    sun          = ephem.Sun()
    emayall.date = time.datetime

    sun.compute(emayall)

    sunra        = sun.ra  * 180. / np.pi
    sundec       = sun.dec * 180. / np.pi
    sunalt       = sun.alt * 180. / np.pi
    
    ##  sunpos   = SkyCoord(ra = sunra * u.degree, dec = sundec * u.degree, frame='icrs').transform_to(AltAz(obstime=Time(mjd, format='mjd'), location=mayall))
    ##  sunalt   = sunpos.alt.degree

    ## 
    isin         = copy.copy(indices)
    
    pos          = SkyCoord(ra = ra[indices] * u.degree, dec = dec[indices] * u.degree, frame='icrs').transform_to(AltAz(obstime=time, location=mayall))
    alt          = pos.alt.degree
    zen          = pos.zen.degree
    airmass      = pos.secz   

    _isin        = list(range(len(alt))) 
    
    ##  az       = np.array([x.az.value  for x in pos]);  airmass = np.array([x.secz.value for x in pos])  
    ##  alt      = np.array([x.alt.value for x in pos])

    islow        = alt < min_altitude
    toremove     = np.argwhere(islow == True)[:,0]
        
    for index in sorted(toremove, reverse=True):
      isin.pop(index)
      _isin.pop(index)
    
    ##  Calculate the local apparent sidereal time in degrees.                                                                                                                                                 
    LST       = LST0 + dLST * (mjd - MJD0)
    hourangle = LST - ra

    moondec, moonra = desisurvey.desiephem.get_object_interpolator(dat._table[i], 'moon', altaz=False)(mjd)

    moonpos         = SkyCoord(ra = moonra * u.degree, dec = moondec * u.degree, frame='icrs').transform_to(AltAz(obstime=time, location=mayall))
    moonalt         = moonpos.alt.degree
     
    if moonalt >= 0.0:
      too_close     = desisurvey.utils.separation_matrix([moonra] * u.deg, [moondec] * u.deg, ra[isin] * u.deg, dec[isin] * u.deg, avoid_bodies['moon'])[0]

      toremove      = np.argwhere(too_close == True)[:,0]

      for index in sorted(toremove, reverse=True):
        isin.pop(index)
        _isin.pop(index)
        
      moons.append([moonra, moondec])
    
    ##  Is this actually bright time?
    moonseps    = desisurvey.utils.separation_matrix([moonra]   * u.deg,   [moondec] * u.deg, ra[isin] * u.deg, dec[isin] * u.deg, None)[0]
    sunseps     = desisurvey.utils.separation_matrix([sunra] * u.deg, [sundec] * u.deg, ra[isin] * u.deg, dec[isin] * u.deg, None)[0]

    ## expfacs  = desisurvey.etc.bright_exposure_factor(moonfrac, moonalt, np.array([x.value for x in moonseps]), sunalt, np.array([x.value for x in sunseps]), airmass[_isin])
    
    sunfacs     = sunmodel(np.array([x.value for x in sunseps]), airmass[_isin], sunalt)
    moonfacs    = moonmodel(np.array([x.value for x in moonseps]), moonalt, moonfrac, airmass[_isin])

    ##  Includes airmass dependence. 
    expfacs     = sunfacs * moonfacs
    
    isbright    = np.all(expfacs >= 2.5)
    
    where       = np.where(expfacs == expfacs.min())[0][0]

    ##  hrs
    eff_time    = dt / expfacs[where]
    
    print('{} \t {:.2f} \t {} \t {:.2f} \t {:.2f} \t {:.2f} \t {:.2f} \t {} \t {:.2f} \t {:.2f} \t {:.2f} \t {:.6f} \t {:.1f} \t {:.2f} \t {:.2f} \t {:.2f} \t {:.2f}'.format(time.datetime, LST, program,\
                                                                                                                                                                              moonfrac, moonalt, sunalt, expfacs.min(),\
                                                                                                                                                                              isbright, moonseps[where].value,\
                                                                                                                                                                              sunseps[where].value,\
                                                                                                                                                                              airmass[_isin][where].value,\
                                                                                                                                                                              eff_time, mjd, moonra, moondec, sunra,\
                                                                                                                                                                              sundec))          
    assert  expfacs.min() > 0.0

    ##
    hrs_visible[i, isin, program] += dt

    brights.append([program, moonfrac, moonalt, sunalt, expfacs.min(), isbright, moonseps[where].value, sunseps[where].value, airmass[_isin][where].value, eff_time, mjd, moonra, moondec, sunra, sundec, LST])

  
  for program in range(3):
    np.savetxt('visibility-{}/visibility-nofullmoon-{}-{}.txt'.format(prefix, nnights, program), hrs_visible[output, :, program])

##
np.savetxt('visibility-{}/moons.txt'.format(prefix), np.array(moons), fmt='%.4lf')

with open('visibility-{}/noons.txt'.format(prefix), 'w') as f:
    for item in isonoons:
        f.write("%s\n" % item)

##
isbright = np.array(isbright)

np.savetxt('brights.txt', brights, fmt='%.4lf')

print('\n\nDone.\n\n')

