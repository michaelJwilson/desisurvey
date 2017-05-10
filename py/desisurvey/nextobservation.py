"""Select the next tile to observe during a night.
"""
from __future__ import print_function, division

import warnings

import numpy as np

import astropy.io.fits
import astropy.time
import astropy.coordinates
import astropy.utils.exceptions
import astropy.units as u

import desiutil.log

import desisurvey.avoidobject
import desisurvey.utils
import desisurvey.config


def nextFieldSelector(obsplan, mjd, progress):
    """Select the next tile to observe during a night.

    Returns the first tile for which the current time falls inside
    its assigned LST window and is far enough from the Moon and
    planets.

    Parameters
    ----------
    obsplan : string
        Name of FITS file containing the afternoon plan
    mjd : float
        Current MJD for selecting the next tile to observe.
    progress : desisurvey.progress.Progress
        Record of observations made so far.  Will not be modified by
        calling this function.

    Returns
    -------
    dict
        Dictionary describing the next tile to observe or None if no
        suitable target is available.  The dictionary will contain the
        following keys: tileID, RA, DEC, Program, Ebmv, moon_illum_frac,
        MoonDist, MoonAlt and overhead.  Overhead is the delay (with time
        units) before the shutter can be opened due to slewing and reading out
        any previous exposure.
    """
    log = desiutil.log.get_logger()

    # Look up configuration parameters.
    config = desisurvey.config.Configuration()
    max_moondist = config.avoid_bodies.moon().to(u.deg).value
    max_zenith_angle = 90 * u.deg - config.min_altitude()

    # Read the afternoon plan.
    hdulist = astropy.io.fits.open(obsplan)
    tiledata = hdulist[1].data
    num_tiles = len(tiledata)
    moonfrac = hdulist[1].header['MOONFRAC']

    # Convert the LST ranges for each tile from hours to degrees.
    tmin = tiledata['LSTMIN'] * 15
    tmax = tiledata['LSTMAX'] * 15

    # Look up the plan exposure time and convert to degrees.
    lst_explen = tiledata['EXPLEN'] / 240.0

    # Convert the current MJD to LST and an astropy time.
    lst = desisurvey.utils.mjd2lst(mjd)
    when = astropy.time.Time(mjd, format='mjd')

    # Initialize pointings for each possible tile.
    ra = tiledata['RA']
    dec = tiledata['DEC']
    proposed = astropy.coordinates.SkyCoord(ra=ra * u.deg, dec=dec * u.deg)

    # Calculate current zenith angle for each possible tile.
    obs = desisurvey.utils.get_observer(when, alt=90 * u.deg, az=0 * u.deg)
    zenith = obs.transform_to(astropy.coordinates.ICRS)
    zenith_angles = proposed.separation(zenith)

    # Determine the previous pointing if we need to include slew time
    # in the overhead calcluations.
    if progress.last_tile is None:
        # No slew needed for next exposure.
        previous = None
        # No readout time needed for previous exposure.
        deadtime = config.readout_time()
    else:
        last = progress.last_tile
        # Where was the telescope last pointing?
        previous = astropy.coordinates.SkyCoord(
            ra=last['ra'] * u.deg, dec=last['dec'] * u.deg)
        # How much time has elapsed since the last exposure ended?
        last_end = (last['mjd'] + last['exptime'] / 86400.).max()
        deadtime = (mjd - last_end) * u.day

    # Calculate the overhead times in seconds for each possible tile.
    overheads = desisurvey.utils.get_overhead_time(previous, proposed, deadtime)
    # Convert to degrees.
    lst_overheads = overheads.to(u.s).value / 240

    found = False
    for i in range(num_tiles):
        tileID = tiledata['TILEID'][i]
        # Estimate this tile's exposure midpoint LST in the range [0,360] deg.
        lst_midpoint = lst + lst_overheads[i] + 0.5 * lst_explen[i]
        if lst_midpoint >= 360:
            lst_midpoint -= 360
        # Skip tiles whose exposure midpoint falls outside their LST window.
        if lst_midpoint < tmin[i] or lst_midpoint > tmax[i]:
            continue
        # Skip a tile that is currently too close to the horizon.
        if zenith_angles[i] > max_zenith_angle:
            log.info(
                'Tile {0} is too close to the horizon ({1:.f})'
                .format(tileID, 90 * u.deg - zenith_angles[i]))
            continue
        # Calculate the moon separation and altitude angles.
        moondist, moonalt, _ = desisurvey.avoidobject.moonLoc(
            when, ra[i], dec[i])
        # Skip a tile that is currently too close to a moon above the horizon.
        if moonalt > 0 and moondist <= max_moondist:
            log.debug(
                'Tile {0} is too close to moon with alt={1:.1f}, sep={2:.1f}.'
                .format(tileID, moonalt, moondist))
            continue
        # Skip a tile that is currently too close to any planets.
        if not desisurvey.avoidobject.avoidObject(when, ra[i], dec[i]):
            log.info('Tile {0} is too close to a planet.'.format(tileID))
            continue
        # Does this tile still needs more exposure?
        if progress.get_tile(tileID)['status'] < 2:
            found = True
            break

    if found:
        # Return a dictionary of parameters needed to observe this tile.
        tile = tiledata[i]
        target = {'tileID' : tile['TILEID'], 'RA' : tile['RA'],
                  'DEC' : tile['DEC'], 'Program': tile['PROGRAM'],
                  'Ebmv' : tile['EBV_MED'], 'moon_illum_frac': moonfrac,
                  'MoonDist': moondist, 'MoonAlt': moonalt,
                  'overhead': overheads[i]}
    else:
        target = None

    return target
