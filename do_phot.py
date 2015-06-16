import os, re
import pyfits as fits
from pylab import *

f_ndx = '2015-03-09-0227'

reg = 'photometry/((\w+)-(\w+)-(\w+)-(\w+)_(\w+)_(\w+)).fits'

def do_phot(f, f_r=None):
    cat_name = 'photometry/' + f.split('/')[1].split('.')[0] + '.cat'
    chk_name = 'photometry/' + f.split('/')[1].split('.')[0] + '_chk.fits'

    if f_r==None:
        cmd = 'sex ' + f + ' -CATALOG_NAME ' + cat_name \
            + ' -c sexf.config'
    else:
        cmd = 'sex ' + f + ' -CATALOG_NAME ' + cat_name \
            + ' -CHECKIMAGE_NAME ' + chk_name + \
            ' -ASSOC_NAME ' + ' photometry/%s_asc.txt ' % f_ndx + ' -c sexf2.config'

    os.system(cmd)

def run_do_phot(patt):
    lc1 = []
    lc2 = []
    tt = []
    files = os.popen('ls ' + patt).read().split()
    for f in files:
        ma = re.match(reg, f)
        if ma.group(7) != '00':
            do_phot(f, 1)
            time = fits.open(f)[0].header['MJD']
            data = loadtxt('photometry/' + f.split('/')[1].split('.')[0] + '.cat')
            if data.shape == (2, 4):
                ndx1 = where(data[:,0] == 1)[0]
                ndx2 = where(data[:,0] == 2)[0]
                if ndx1.shape[0]==1 & ndx2.shape[0]==1:
                    lc1.append(data[ndx1,3])
                    lc2.append(data[ndx2,3])
                    tt.append(time)
                print size(lc1), size(lc2)

    ff = open('lc_%s.txt' % f_ndx , 'w')

    for i in range(size(tt)):
        ff.write('%20.8f %7.4f %7.4f\n' % (tt[i], lc1[i],lc2[i]))

    ff.close()


run_do_phot('photometry/%s*.fits' % f_ndx)
