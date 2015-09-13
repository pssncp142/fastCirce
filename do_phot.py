import os, re
import pyfits as fits
from pylab import *

f_ndx = '2015-06-17-0064'

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
            ' -ASSOC_NAME ' + ' asc/%s_asc.txt ' % f_ndx + ' -c sexf2.config'

    os.system(cmd)

def run_do_phot(patt):
    asc = len(open('asc/%s_asc.txt' % f_ndx).read().split('\n')[:-1])
    print open('asc/%s_asc.txt' % f_ndx).read().split('\n')
    lcs = []
    for i in arange(asc): lcs.append([])
    print lcs
    tt = []
    files = os.popen('ls ' + patt).read().split()
    for f in files:
        ma = re.match(reg, f)
        if ma.group(7) != '00':
            do_phot(f, 1)
            time = fits.open(f)[0].header['MJD']
            data = loadtxt('photometry/' + f.split('/')[1].split('.')[0] + '.cat')
            for i in arange(asc):
                ndx = where(data[:,0] == i+1)[0]
                if ndx.shape[0] == 1:
                    lcs[i].append(data[ndx,3])
                else:
                    lcs[i].append(-8)
            tt.append(time)

    ff = open('lc_%s.txt' % f_ndx , 'w')

    print lcs
    print tt

    for i in range(size(tt)):
        ff.write('%.8f' % tt[i])
        for j in range(asc):
            print j, i
            ff.write(',%.4f ' % (lcs[j][i]))
        ff.write('\n')

    ff.close()


run_do_phot('photometry/%s*.fits' % f_ndx)
