import numpy         as     np

from   astropy.io    import fits
from   astropy.table import Table, Column, vstack, join


def expderived(exposures, tiles=None, _mjd='MJD', unmasked=True):
    '''                                                                                                                                                                                                                                                                                                                                                                                                                         
    Given RA, DEC, and MJD from exposures table, return derived proprs for Kitt Peak, e.g. MOONFRAC, MOONSEP etc.                                                                                                                                                                                                                                                                                                               
    Note:  temporary copy of https://github.com/changhoonhahn/feasiBGS/blob/master/run/bright_exposure/surveysim_output.py, line 384.                                                                                                                                                                                                                                                                                           
    '''
    import ephem
    import desisurvey.config
    import desisurvey.utils   as      dutils
    import astropy.units      as      u

    from   astropy.time       import  Time


    config                = desisurvey.config.Configuration()

    mayall                = ephem.Observer()
    mayall.lat            = config.location.latitude().to(u.rad).value
    mayall.lon            = config.location.longitude().to(u.rad).value
    mayall.elevation      = config.location.elevation().to(u.m).value

    ##  Observed time (MJD).                                                                                                                                                                                                                                                                                                                                                                                                    
    mjd                   = Time(exposures[_mjd].quantity.value, format='mjd')

    exposures['MOONALT']  = np.zeros(len(mjd))
    exposures['MOONRA']   = np.zeros(len(mjd))
    exposures['MOONDEC']  = np.zeros(len(mjd))
    exposures['MOONFRAC'] = np.zeros(len(mjd))

    exposures['SUNALT']   = np.zeros(len(mjd))
    exposures['SUNRA']    = np.zeros(len(mjd))
    exposures['SUNDEC']   = np.zeros(len(mjd))

    _moon = ephem.Moon()
    _sun  = ephem.Sun()

    print('Solving for lunar positions.')
    
    for i in range(len(mjd)):
        mayall.date              = mjd.datetime[i]

        _moon.compute(mayall)
        _sun.compute(mayall)

        exposures['MOONALT'][i]  = _moon.alt
        exposures['MOONRA'][i]   = _moon.ra
        exposures['MOONDEC'][i]  = _moon.dec
        exposures['MOONFRAC'][i] = _moon.moon_phase

        exposures['SUNALT'][i]   = _sun.alt
        exposures['SUNRA'][i]    = _sun.ra
        exposures['SUNDEC'][i]   = _sun.dec

    ##  
    exposures['MOONALT'] *= 180. / np.pi
    exposures['MOONRA']  *= 180. / np.pi
    exposures['MOONDEC'] *= 180. / np.pi

    exposures['SUNALT']  *= 180. / np.pi
    exposures['SUNRA']   *= 180. / np.pi
    exposures['SUNDEC']  *= 180. / np.pi

    if tiles is not None:
        print('\n\nJoining with Tiles file.')

        exposures            = join(exposures, tiles, keys=['TILEID'], join_type='outer', table_names=['', '~'], uniq_col_name='{col_name}{table_name}')
        exposures.remove_column('AIRMASS~')

        exposures['MOONSEP'] = np.diag(dutils.separation_matrix(exposures['MOONRA'].quantity * u.deg, exposures['MOONDEC'].quantity * u.deg,\
                                       np.atleast_1d(exposures['RA'].quantity * u.deg), np.atleast_1d(exposures['DEC'].quantity * u.deg)))

        exposures['SUNSEP']  = np.diag(dutils.separation_matrix(exposures['SUNRA'].quantity * u.deg, exposures['SUNDEC'].quantity * u.deg,\
                                       np.atleast_1d(exposures['RA'].quantity * u.deg), np.atleast_1d(exposures['DEC'].quantity * u.deg)))

    if unmasked:
      exposures.remove_rows(np.where(exposures['EXPTIME'].mask == True))

    return  exposures
