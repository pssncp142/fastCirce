import pyfits as fits 
import pywcs

f_name = '/scratch1/ydallilar/2015-06-17/CIRCE2015-06-17-0082.fits'

hdul = fits.open(f_name)

head = hdul[0].header

wcs = pywcs.WCS(head)

nhead = wcs.to_header()

nhead['CRPIX2'] = head['CRPIX2'] - 0*head['W_Y_BEG']

hdu = fits.PrimaryHDU(hdul[10].data,header=nhead)
hdu.writeto('test1-4.fits', clobber=True)
