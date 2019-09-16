import astropy.io.fits as fits

from   desimodel import io
from   astropy.table import Table


##  tfile         = io.findfile('footprint/desi-tiles.fits')
tfile             = '/global/homes/m/mjwilson/repos/chh/desisurvey/py/desisurvey/desi-tiles.fits'

tiles             = fits.open(tfile)

tiles             = Table(tiles[1].data)
tiles             = tiles[tiles['PROGRAM'] == 'BRIGHT']
tiles             = tiles[tiles['IN_DESI']  >        0]

tiles.sort('TILEID')

##  Have bright passes start at zero. 
tiles['PASS']    -= 5

print(tiles)

tiles.write('bright-tiles.fits', format='fits', overwrite=True)
