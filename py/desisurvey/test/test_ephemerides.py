import unittest
import os
import uuid
import datetime

import numpy as np

from astropy.time import Time
from astropy.coordinates import ICRS, AltAz
import astropy.units as u
import astropy.io

from desisurvey.ephemerides import Ephemerides, get_grid, get_object_interpolator
from desisurvey.config import Configuration
from desisurvey.utils import get_date, get_location


class TestEphemerides(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.origdir = os.getcwd()
        cls.testdir = os.path.abspath('./test-{}'.format(uuid.uuid4()))
        os.mkdir(cls.testdir)
        os.chdir(cls.testdir)
        # Configure a CSV reader for the Horizons output format.
        csv_reader = astropy.io.ascii.Csv()
        csv_reader.header.comment = r'[^ ]'
        csv_reader.data.start_line = 35
        csv_reader.data.end_line = 203
        # Read moon ephemerides for the first week of 2020.
        path = astropy.utils.data._find_pkg_data_path(
            os.path.join('data', 'horizons_2020_week1_moon.csv'),
            package='desisurvey')
        cls.table = csv_reader.read(path)
        # Horizons CSV file has a trailing comma on each line.
        cls.table.remove_column('col10')
        # Use more convenient column names.
        names = ('date', 'jd', 'sun', 'moon', 'ra', 'dec',
                 'az', 'alt', 'lst', 'frac')
        for old_name, new_name in zip(cls.table.colnames, names):
            cls.table[old_name].name = new_name

    @classmethod
    def tearDownClass(cls):
        os.chdir(cls.origdir)
        if os.path.isdir(cls.testdir):
            import shutil
            shutil.rmtree(cls.testdir)

    def test_getephem(self):
        """Tabulate one month of ephemerides"""
        start = datetime.date(2019, 9, 1)
        stop = datetime.date(2019, 10, 1)
        ephem = Ephemerides(start, stop, use_cache=False, write_cache=False)
        self.assertEqual(ephem.start.mjd, ephem.get_row(0)['noon'])
        self.assertEqual(ephem.start.mjd, ephem.get_night(start)['noon'])
        self.assertEqual(ephem.stop.mjd, ephem.get_row(-1)['noon'] + 1)
        self.assertEqual(ephem.start.mjd,
                         ephem.get_night(ephem.start)['noon'])
        self.assertEqual(ephem.num_nights,
                         int(round(ephem.stop.mjd - ephem.start.mjd)))

        etable = ephem._table
        self.assertEqual(len(etable), 30)
        self.assertTrue(np.all(etable['dusk'] > etable['noon']))
        self.assertTrue(np.all(etable['dawn'] > etable['dusk']))
        self.assertTrue(np.all(etable['dusk'] > etable['brightdusk']))
        self.assertTrue(np.all(etable['dawn'] < etable['brightdawn']))
        self.assertGreater(np.max(etable['moon_illum_frac']), 0.99)
        self.assertLessEqual(np.max(etable['moon_illum_frac']), 1.0)
        self.assertLess(np.min(etable['moon_illum_frac']), 0.01)
        self.assertGreaterEqual(np.min(etable['moon_illum_frac']), 0.00)
        self.assertTrue(np.all(etable['moonrise'] < etable['moonset']))

        for i in range(ephem.num_nights):

            x = ephem.get_row(i)
            date = Time(x['noon'], format='mjd').datetime.date()
            night = date.strftime('%Y%m%d')
            for key in [
                    'brightdusk', 'brightdawn',
                    'dusk', 'dawn',
                ]:
                #- AZ local time
                localtime = Time(x[key], format='mjd') - 7*u.hour
                #- YEARMMDD of sunset for that time
                yearmmdd = (localtime - 12*u.hour).to_datetime().strftime('%Y%m%d')
                msg = '{} != {} for {}={}'.format(night, yearmmdd, key, x[key])
                self.assertEqual(night, yearmmdd, msg)

    def test_full_moon(self):
        """Verify that the full moon break in Sep-2019 occurs on days 10-16"""
        start = datetime.date(2019, 9, 1)
        stop = datetime.date(2019, 9, 30)
        ephem = Ephemerides(start, stop, use_cache=False, write_cache=False)
        full = np.empty(ephem.num_nights, bool)
        for i in range(ephem.num_nights):
            night = start + datetime.timedelta(days=i)
            full[i] = ephem.is_full_moon(night)
        expected = np.zeros_like(full, bool)
        expected[9:16] = True
        self.assertTrue(np.all(full == expected))

    def test_get_grid(self):
        """Verify grid calculations"""
        for step_size in (1 * u.min, 0.3 * u.hour):
            for night_start in (-6 * u.hour, -6.4 * u.hour):
                g = get_grid(step_size, night_start)
                self.assertTrue(g[0] == night_start.to(u.day).value)
                self.assertAlmostEqual(g[1] - g[0], step_size.to(u.day).value)
                self.assertAlmostEqual(g[-1] - g[0],
                                (len(g) - 1) * step_size.to(u.day).value)

    def test_moon_phase(self):
        """Verfify moon illuminated fraction for first week of 2020"""
        ephem = Ephemerides(
            get_date('2019-12-31'), get_date('2020-02-02'),
            use_cache=False, write_cache=False)
        for i, jd in enumerate(self.table['jd']):
            t = Time(jd, format='jd')
            frac = ephem.get_moon_illuminated_fraction(t.mjd)
            truth = 1e-2 * self.table['frac'][i]
            self.assertTrue(abs(frac - truth) < 0.01)

    def test_moon_radec(self):
        """Verify moon (ra,dec) for first week of 2020"""
        ephem = Ephemerides(
            get_date('2019-12-31'), get_date('2020-02-02'),
            use_cache=False, write_cache=False)
        for i, jd in enumerate(self.table['jd']):
            t = Time(jd, format='jd')
            night = ephem.get_night(t)
            f_moon = get_object_interpolator(night, 'moon', altaz=False)
            dec, ra = f_moon(t.mjd)
            truth = ICRS(ra=self.table['ra'][i] * u.deg,
                         dec=self.table['dec'][i] * u.deg)
            calc = ICRS(ra=ra * u.deg, dec=dec * u.deg)
            sep = truth.separation(calc)
            self.assertTrue(abs(sep.to(u.deg).value) < 0.3)

    def test_moon_altaz(self):
        """Verify moon (alt,az) for first week of 2020"""
        ephem = Ephemerides(
            get_date('2019-12-31'), get_date('2020-02-02'),
            use_cache=False, write_cache=False)
        location = get_location()
        for i, jd in enumerate(self.table['jd']):
            t = Time(jd, format='jd')
            night = ephem.get_night(t)
            f_moon = get_object_interpolator(night, 'moon', altaz=True)
            alt, az = f_moon(t.mjd)
            truth = AltAz(alt=self.table['alt'][i] * u.deg,
                          az=self.table['az'][i] * u.deg,
                          obstime=t, location=location, pressure=0)
            calc = AltAz(alt=alt * u.deg, az=az * u.deg,
                         obstime=t, location=location, pressure=0)
            sep = truth.separation(calc)
            self.assertTrue(abs(sep.to(u.deg).value) < 0.3)


if __name__ == '__main__':
    unittest.main()
