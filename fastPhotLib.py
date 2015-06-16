
import pyfits as fits
from numpy import *
import os, re

frame_mod = '3-1'
path = 'reduced_data/'
path_out = 'photometry/'

f_name = ['CIRCE2015-03-09-0227.fits']

patt = 'CIRCE((\w+)-(\w+)-(\w+)-(\w+)).fits'

def prepare_phot(name):

    if not os.path.isdir(path_out):
        os.mkdir(path_out)
    
    ma = re.match(patt, name)
    f_ndx = ma.group(1)
    
    ff = os.popen('ls %smod%s_%s_*.fits' % (path, frame_mod, f_ndx)).\
         read().split()

    for i,f in enumerate(ff):
        print i, f
        hdul = fits.open(f)
        ndx = 1
        phdu = fits.PrimaryHDU(hdul[-1].data,header=hdul[-1].header)
        phdu.writeto('%s%s_%02d_%02d.fits' % \
                     (path_out, f_ndx, i+1, 0),
                     clobber=True)

        for hdu in hdul[1:-1]:
            phdu = fits.PrimaryHDU(hdu.data,header=hdu.header)
            phdu.writeto('%s%s_%02d_%02d.fits' % \
                         (path_out, f_ndx, i+1, ndx),
                         clobber=True)
            ndx += 1
            

            
if __name__ == '__main__':
    prepare_phot(f_name[0])
