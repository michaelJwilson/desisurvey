"""Utility functions for survey planning and scheduling.
"""
from __future__ import print_function, division

import datetime
import os
import warnings

import numpy as np

import pytz

import astropy.time
import astropy.coordinates
import astropy.utils.iers
import astropy.utils.data
import astropy.utils.exceptions
import astropy._erfa.core
import astropy.units as u

import desiutil.log

import desisurvey.config


_telescope_location = None
_iers_is_frozen = False
_dome_closed_probabilities = None


def dome_closed_probabilities():
    """Return an array of monthly dome-closed probabilities.

    Returns
    -------
    array
        Array of 12 probabilities in the range 0-1.
    """
    global _dome_closed_probabilities
    if _dome_closed_probabilities is None:
        config = desisurvey.config.Configuration()
        _dome_closed_probabilities = np.empty(12)
        for i, month in enumerate(('jan', 'feb', 'mar', 'apr', 'may', 'jun',
                                   'jul', 'aug', 'sep', 'oct', 'nov', 'dec')):
            _dome_closed_probabilities[i] = (
                # Convert from percentage to fraction.
                getattr(config.dome_closed_probability, month)() / 100.)
    return _dome_closed_probabilities


def freeze_iers(name='iers_frozen.ecsv', ignore_warnings=True):
    """Use a frozen IERS table saved with this package.

    This should be called at the beginning of a script that calls
    astropy time and coordinates functions which refer to the UT1-UTC
    and polar motions tabulated by IERS.  The purpose is to ensure
    identical results across systems and astropy releases, to avoid a
    potential network download, and to eliminate some astropy warnings.

    After this call, the loaded table will be returned by
    :func:`astropy.utils.iers.IERS_Auto.open()` and treated like a
    a normal IERS table by all astropy code.  Specifically, this method
    registers an instance of a custom IERS_Frozen class that inherits from
    IERS_B and overrides
    :meth:`astropy.utils.iers.IERS._check_interpolate_indices` to prevent
    any IERSRangeError being raised.

    See http://docs.astropy.org/en/stable/utils/iers.html for details.

    This function returns immediately after the first time it is called,
    so it it safe to insert anywhere that consistent IERS models are
    required, and subsequent calls with different args will have no
    effect.

    The :func:`desisurvey.utils.plot_iers` function is useful for inspecting
    IERS tables and how they are extrapolated to DESI survey dates.

    Parameters
    ----------
    name : str
        Name of the file to load the frozen IERS table from. Should normally
        be relative and then refers to this package's data/ directory.
        Must end with the .ecsv extension.
    ignore_warnings : bool
        Ignore ERFA and IERS warnings about future dates generated by
        astropy time and coordinates functions. Specifically, ERFA warnings
        containing the string "dubious year" are filtered out, as well
        as AstropyWarnings related to IERS table extrapolation.
    """
    log = desiutil.log.get_logger()
    if desisurvey.utils._iers_is_frozen:
        log.debug('IERS table already frozen.')
        return
    log.info('Freezing IERS table used by astropy time, coordinates.')

    # Validate the save_name extension.
    _, ext = os.path.splitext(name)
    if ext != '.ecsv':
        raise ValueError('Expected .ecsv extension for {0}.'.format(name))

    # Locate the file in our package data/ directory.
    if not os.path.isabs(name):
        name = astropy.utils.data._find_pkg_data_path(
            os.path.join('data', name))
    if not os.path.exists(name):
        raise ValueError('No such IERS file: {0}.'.format(name))

    # Clear any current IERS table.
    astropy.utils.iers.IERS.close()
    # Initialize the global IERS table. We load the table by
    # hand since the IERS open() method hardcodes format='cds'.
    try:
        table = astropy.table.Table.read(name, format='ascii.ecsv').filled()
    except IOError:
        raise RuntimeError('Unable to load IERS table from {0}.'.format(name))

    # Define a subclass of IERS_B that overrides _check_interpolate_indices
    # to prevent any IERSRangeError being raised.
    class IERS_Frozen(astropy.utils.iers.IERS_B):
        def _check_interpolate_indices(self, indices_orig, indices_clipped,
                                       max_input_mjd): pass

    # Create and register an instance of this class from the table.
    iers = IERS_Frozen(table)
    astropy.utils.iers.IERS.iers_table = iers
    # Prevent any attempts to automatically download updated IERS-A tables.
    astropy.utils.iers.conf.auto_download = False
    astropy.utils.iers.conf.auto_max_age = None
    astropy.utils.iers.conf.iers_auto_url = 'frozen'
    # Sanity check.
    if not (astropy.utils.iers.IERS_Auto.open() is iers):
        raise RuntimeError('Frozen IERS is not installed as the default.')

    if ignore_warnings:
        warnings.filterwarnings(
            'ignore', category=astropy._erfa.core.ErfaWarning, message=
            r'ERFA function \"[a-z0-9_]+\" yielded [0-9]+ of \"dubious year')
        warnings.filterwarnings(
            'ignore', category=astropy.utils.exceptions.AstropyWarning,
            message=r'Tried to get polar motions for times after IERS data')
        warnings.filterwarnings(
            'ignore', category=astropy.utils.exceptions.AstropyWarning,
            message=r'\(some\) times are outside of range covered by IERS')

    # Shortcircuit any subsequent calls to this function.
    desisurvey.utils._iers_is_frozen = True


def update_iers(save_name='iers_frozen.ecsv', num_avg=1000):
    """Update the IERS table used by astropy time, coordinates.

    Downloads the current IERS-A table, replaces the last entry (which is
    repeated for future times) with the average of the last ``num_avg``
    entries, and saves the table in ECSV format.

    This should only be called every few months, e.g., with major releases.
    The saved file should then be copied to this package's data/ directory
    and committed to the git repository.

    Requires a network connection in order to download the current IERS-A table.
    Prints information about the update process.

    The :func:`desisurvey.utils.plot_iers` function is useful for inspecting
    IERS tables and how they are extrapolated to DESI survey dates.

    Parameters
    ----------
    save_name : str
        Name where frozen IERS table should be saved. Must end with the
        .ecsv extension.
    num_avg : int
        Number of rows from the end of the current table to average and
        use for calculating UT1-UTC offsets and polar motion at times
        beyond the table.
    """
    # Validate the save_name extension.
    _, ext = os.path.splitext(save_name)
    if ext != '.ecsv':
        raise ValueError('Expected .ecsv extension for {0}.'.format(save_name))

    # Download the latest IERS_A table
    iers = astropy.utils.iers.IERS_A.open(astropy.utils.iers.IERS_A_URL)
    last = astropy.time.Time(iers['MJD'][-1], format='mjd').datetime
    print('Updating to current IERS-A table with coverage up to {0}.'
          .format(last.date()))

    # Loop over the columns used by the astropy IERS routines.
    for name in 'UT1_UTC', 'PM_x', 'PM_y':
        # Replace the last entry with the mean of recent samples.
        mean_value = np.mean(iers[name][-num_avg:].value)
        unit = iers[name].unit
        iers[name][-1] = mean_value * unit
        print('Future {0:7s} = {1:.3}'.format(name, mean_value * unit))

    # Strip the original table metadata since ECSV cannot handle it.
    # We only need a single keyword that is checked by IERS_Auto.open().
    iers.meta = dict(data_url='frozen')

    # Save the table. The IERS-B table provided with astropy uses the
    # ascii.cds format but astropy cannot write this format.
    iers.write(save_name, format='ascii.ecsv')
    print('Wrote updated table to {0}.'.format(save_name))


def get_overhead_time(current_pointing, new_pointing, deadtime=0):
    """Compute the instrument overhead time between exposures.

    Use a model of the time required to slew and focus, in parallel with
    reading out the previous exposure.

    With no slew or readout required, the minimum overhead is set by the
    time needed to focus the new exposure.

    The calculation will be automatically broadcast over an array of new
    pointings and return an array of overhead times.

    Parameters
    ----------
    current_pointing : :class:`astropy.coordinates.SkyCoord` or None
        Current pointing of the telescope.  Do not include any slew overhead
        when None.
    new_pointing : :class:`astropy.coordinates.SkyCoord`
        New pointing(s) of the telescope.
    deadtime : :class:`astropy.units.Quantity`, optional
        Amount of deadtime elapsed since end of any previous exposure.
        Used to ensure that the overhead time is sufficient to finish
        reading out the previous exposure. Must be >= 0.
        Defaults to zero seconds.

    Returns
    -------
    :class:`astropy.units.Quantity`
        Overhead time(s) for each new_pointing.
    """
    if not isinstance(deadtime, u.Quantity):
        deadtime = deadtime * u.s
    if deadtime.to(u.s).value < 0:
        raise ValueError('Expected deadtime >= 0 (got {0}).'.format(deadtime))
    config = desisurvey.config.Configuration()
    if current_pointing is not None:
        # Calculate the amount that each axis needs to move in degrees.
        # The ra,dec attributes of a SkyCoord are always in the ranges
        # [0,360] and [-90,+90] degrees.
        delta_dec = np.fabs(
            (new_pointing.dec - current_pointing.dec).to(u.deg).value)
        delta_ra = np.fabs(
            (new_pointing.ra - current_pointing.ra).to(u.deg).value)
        # Handle wrap around in RA
        delta_ra = 180 - np.fabs(delta_ra - 180)
        # The slew time is determined by the axis motor with the most travel.
        max_travel = np.maximum(delta_ra, delta_dec) * u.deg
        moving = max_travel > 0
        overhead = max_travel / config.slew_rate()
        overhead[moving] += config.slew_overhead()
    else:
        overhead = np.zeros(new_pointing.shape) * u.s
    # Add the constant focus time.
    overhead += config.focus_time()
    # Overhead is at least the remaining readout time for the last exposure.
    overhead = np.maximum(overhead, config.readout_time() - deadtime)

    return overhead


def get_location():
    """Return the telescope's earth location.

    The location object is cached after the first call, so there is no need
    to cache this function's return value externally.

    Returns
    -------
    astropy.coordinates.EarthLocation
    """
    global _telescope_location
    if _telescope_location is None:
        config = desisurvey.config.Configuration()
        _telescope_location = astropy.coordinates.EarthLocation.from_geodetic(
            lat=config.location.latitude(),
            lon=config.location.longitude(),
            height=config.location.elevation())
    return _telescope_location


def get_observer(when, alt=None, az=None):
    """Return the AltAz frame for the telescope at the specified time(s).

    Refraction corrections are not applied (for now).

    The returned object is automatically broadcast over input arrays.

    Parameters
    ----------
    when : astropy.time.Time
        One or more times when the AltAz transformations should be calculated.
    alt : astropy.units.Quantity or None
        Local altitude angle(s)
    az : astropy.units.Quantity or None
        Local azimuth angle(s)

    Returns
    -------
    astropy.coordinates.AltAz
        AltAz frame object suitable for transforming to/from local horizon
        (alt, az) coordinates.
    """
    if alt is not None and az is not None:
        kwargs = dict(alt=alt, az=az)
    elif alt is not None or az is not None:
        raise ValueError('Must specify both alt and az.')
    else:
        kwargs = {}
    return astropy.coordinates.AltAz(
        location=get_location(), obstime=when, pressure=0, **kwargs)


def cos_zenith_to_airmass(cosZ):
    """Convert a zenith angle to an airmass.

    Uses the Rozenberg 1966 interpolation formula, which gives reasonable
    results for high zenith angles, with a horizon air mass of 40.
    https://en.wikipedia.org/wiki/Air_mass_(astronomy)#Interpolative_formulas
    Rozenberg, G. V. 1966. "Twilight: A Study in Atmospheric Optics."
    New York: Plenum Press, 160.

    The value of cosZ is clipped to [0,1], so observations below the horizon
    return the horizon value (~40).

    Parameters
    ----------
    cosZ : float or array
        Cosine of angle(s) to convert.

    Returns
    -------
    float or array
        Airmass value(s) >= 1.
    """
    cosZ = np.clip(np.asarray(cosZ), 0., 1.)
    return np.clip(1. / (cosZ + 0.025 * np.exp(-11 * cosZ)), 1., None)


def get_airmass(when, ra, dec):
    """Return the airmass of (ra,dec) at the specified observing time.

    Uses :func:`cos_zenith_to_airmass`.

    Parameters
    ----------
    when : astropy.time.Time
        Observation time, which specifies the local zenith.
    ra : astropy.units.Quantity
        Target RA angle(s)
    dec : astropy.units.Quantity
        Target DEC angle(s)
    Returns
    -------
    array or float
        Value of the airmass for each input (ra,dec).
    """
    target = astropy.coordinates.ICRS(ra=ra, dec=dec)
    zenith = get_observer(when, alt=90 * u.deg, az=0 * u.deg
                          ).transform_to(astropy.coordinates.ICRS)
    # Calculate zenith angle in degrees.
    zenith_angle = target.separation(zenith)
    # Convert to airmass.
    return cos_zenith_to_airmass(np.cos(zenith_angle))


def cos_zenith(ha, dec, latitude=None):
    """Calculate cos(zenith) for specified hour angle, DEC and latitude.

    Combine with :func:`cos_zenith_to_airmass` to calculate airmass.

    Parameters
    ----------
    ha : astropy.units.Quantity
        Hour angle(s) to use, with units convertible to angle.
    dec : astropy.units.Quantity
        Declination angle(s) to use, with units convertible to angle.
    latitude : astropy.units.Quantity or None
        Latitude angle to use, with units convertible to angle.
        Defaults to the latitude of :func:`get_location` if None.

    Returns
    -------
    numpy array
        cosine of zenith angle(s) corresponding to the inputs.
    """
    if latitude is None:
        # Use the observatory latitude by default.
        latitude = desisurvey.config.Configuration().location.latitude()
    # Calculate sin(altitude) = cos(zenith).
    cosZ = (np.sin(dec) * np.sin(latitude) +
            np.cos(dec) * np.cos(latitude) * np.cos(ha))
    # Return a plain array (instead of a unitless Quantity).
    return cosZ.value


def is_monsoon(night):
    """Test if this night's observing falls in the monsoon shutdown.

    Uses the monsoon date range defined in the
    :class:`desisurvey.config.Configuration`.  Based on (month, day) comparisons
    rather than day-of-year comparisons, so the monsoon always starts on the
    same calendar date, even in leap years.

    Parameters
    ----------
    night : date
        Converted to a date using :func:`desisurvey.utils.get_date`.

    Returns
    -------
    bool
        True if this night's observing falls during the monsoon shutdown.
    """
    date = get_date(night)

    # Fetch our configuration.
    config = desisurvey.config.Configuration()
    start = config.monsoon_start()
    stop = config.monsoon_stop()

    # Not in monsoon if (day < start) or (day >= stop)
    m, d = date.month, date.day
    if m < start.month or (m == start.month and d < start.day):
        return False
    if m > stop.month or (m == stop.month and d >= stop.day):
        return False

    return True


def local_noon_on_date(day):
    """Convert a date to an astropy time at local noon.

    Local noon is used as the separator between observing nights. The purpose
    of this function is to standardize the boundary between observing nights
    and the mapping of dates to times.

    Generates astropy ErfaWarnings for times in the future.

    Parameters
    ----------
    day : datetime.date
        The day to use for generating a time object.

    Returns
    -------
    astropy.time.Time
        A Time object with the input date and a time corresponding to
        local noon at the telescope.
    """
    # Fetch our configuration.
    config = desisurvey.config.Configuration()

    # Build a datetime object at local noon.
    tz = pytz.timezone(config.location.timezone())
    local_noon = tz.localize(
        datetime.datetime.combine(day, datetime.time(hour=12)))

    # Convert to UTC.
    utc_noon = local_noon.astimezone(pytz.utc)

    # Return a corresponding astropy Time.
    return astropy.time.Time(utc_noon)


def get_date(date):
    """Convert different date specifications into a datetime.date object.

    We use strptime() to convert an input string, so leading zeros are not
    required for strings in the format YYYY-MM-DD, e.g. 2019-8-3 is considered
    valid.

    Instead of testing the input type, we try different conversion methods:
    ``.datetime.date()`` for an astropy time and ``datetime.date()`` for a
    datetime.

    Date specifications that include a time of day (datetime, astropy time, MJD)
    are rounded down to the previous local noon before converting to a date.
    This ensures that all times during a local observing night are mapped to
    the same date, when the night started.  A "naive" (un-localized) datetime
    is assumed to refer to UTC.

    Generates astropy ERFA warnings for future dates.

    Parameters
    ----------
    date : astropy.time.Time, datetime.date, datetime.datetime, string or number
        Specification of the date to return.  A string must have the format
        YYYY-MM-DD (but leading zeros on MM and DD are optional).  A number
        will be interpreted as a UTC MJD value.

    Returns
    -------
    datetime.date
    """
    input_date = date
    # valid types: string, number, Time, datetime, date
    try:
        # Convert bytes to str.
        date = date.decode()
    except AttributeError:
        pass
    try:
        # Convert a string of the form YYYY-MM-DD into a date.
        # This will raise a ValueError for a badly formatted string
        # or invalid date such as 2019-13-01.
        try:
            date = datetime.datetime.strptime(date, '%Y-%m-%d').date()
        except ValueError:
            try:
                date = datetime.datetime.strptime(date, '%Y%m%d').date()
            except ValueError:
                raise
    except TypeError:
        pass
    # valid types: number, Time, datetime, date
    try:
        # Convert a number to an astropy time, assuming it is a UTC MJD value.
        date = astropy.time.Time(date, format='mjd')
    except ValueError:
        pass
    # valid types: Time, datetime, date
    try:
        # Convert an astropy time into a datetime
        date = date.datetime
    except AttributeError:
        pass
    # valid types: datetime, date
    try:
        # Localize a naive datetime assuming it refers to UTC.
        date = pytz.utc.localize(date)
    except (AttributeError, ValueError):
        pass
    # valid types: localized datetime, date
    try:
        # Convert a localized datetime into the date of the previous noon.
        local_tz = pytz.timezone(
            desisurvey.config.Configuration().location.timezone())
        local_time = date.astimezone(local_tz)
        date = local_time.date()
        if local_time.hour < 12:
            date -= datetime.timedelta(days=1)
    except AttributeError:
        pass
    # valid types: date
    if not isinstance(date, datetime.date):
        raise ValueError('Invalid date specification: {0}.'.format(input_date))
    return date


def day_number(date):
    """Return the number of elapsed days since the start of the survey.

    Does not perform any range check that the date is within the nominal
    survey schedule.

    Parameters
    ----------
    date : astropy.time.Time, datetime.date, datetime.datetime, string or number
        Converted to a date using :func:`get_date`.

    Returns
    -------
    int
        Number of elapsed days since the start of the survey.
    """
    config = desisurvey.config.Configuration()
    return (get_date(date) - config.first_day()).days


def separation_matrix(ra1, dec1, ra2, dec2, max_separation=None):
    """Build a matrix of pair-wise separation between (ra,dec) pointings.

    The ra1 and dec1 arrays must have the same shape. The ra2 and dec2 arrays
    must also have the same shape, but it can be different from the (ra1,dec1)
    shape, resulting in a non-square return matrix.

    Uses the Haversine formula for better accuracy at low separations. See
    https://en.wikipedia.org/wiki/Haversine_formula for details.

    Equivalent to using the separations() method of astropy.coordinates.ICRS,
    but faster since it bypasses any units.

    Parameters
    ----------
    ra1 : array
        1D array of n1 RA coordinates in degrees (without units attached).
    dec1 : array
        1D array of n1 DEC coordinates in degrees (without units attached).
    ra2 : array
        1D array of n2 RA coordinates in degrees (without units attached).
    dec2 : array
        1D array of n2 DEC coordinates in degrees (without units attached).
    max_separation : float or None
        When present, the matrix elements are replaced with booleans given
        by (value <= max_separation), which saves some computation.

    Returns
    -------
    array
        Array with shape (n1,n2) with element [i1,i2] giving the 3D separation
        angle between (ra1[i1],dec1[i1]) and (ra2[i2],dec2[i2]) in degrees
        or, if max_separation is not None, booleans (value <= max_separation).
    """
    ra1, ra2 = np.deg2rad(ra1), np.deg2rad(ra2)
    dec1, dec2 = np.deg2rad(dec1), np.deg2rad(dec2)
    if ra1.shape != dec1.shape:
        raise ValueError('Arrays ra1, dec1 must have the same shape.')
    if len(ra1.shape) != 1:
        raise ValueError('Arrays ra1, dec1 must be 1D.')
    if ra2.shape != dec2.shape:
        raise ValueError('Arrays ra2, dec2 must have the same shape.')
    if len(ra2.shape) != 1:
        raise ValueError('Arrays ra2, dec2 must be 1D.')
    havRA12 = 0.5 * (1 - np.cos(ra2 - ra1[:, np.newaxis]))
    havDEC12 = 0.5 * (1 - np.cos(dec2 - dec1[:, np.newaxis]))
    havPHI = havDEC12 + np.cos(dec1)[:, np.newaxis] * np.cos(dec2) * havRA12
    if max_separation is not None:
        # Replace n1 x n2 arccos calls with a single sin call.
        threshold = np.sin(0.5 * np.deg2rad(max_separation)) ** 2
        return havPHI <= threshold
    else:
        return np.rad2deg(np.arccos(np.clip(1 - 2 * havPHI, -1, +1)))
