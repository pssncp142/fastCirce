#------------------------------------------------------------------------------#
# Yigit Dallilar, UF, 05/30/2015                                               #
# Contact: ydallilar@ufl.edu                                                   #
#------------------------------------------------------------------------------#
# CIRCE fast photometry pipeline tools                                         #
# -fastCirceLib.py   version=1.0.1                                             #
# Contains image processing functions                                          #
#------------------------------------------------------------------------------#

import pyfits as fits
import re
import os 
from numpy import *

sep = '-----------------------------------------------------------------------'
patt = 'CIRCE((\w+)-(\w+)-(\w+)-(\w+)).fits'
t_c = 24*60*60

#------------------------------------------------------------------------------#
# Do linearity                                                                 #
#------------------------------------------------------------------------------#

def do_linearity(obj):

    coef = [0.9701104, 1.106431e-5, -7.398100e-10, 2.490898e-14]
    print sep
    print 'Starting linearity correction...'
    dark = fits.open(obj.mst_dark)[0].data
    for f in obj.files:
        hdul = fits.open(obj.path + f)
        print 'Processing ', obj.path + f
        hdul_out = [fits.PrimaryHDU()]
        hdul_out[0].header.update('W_Y_BEG', hdul[0].header['W_Y_BEG'])
        hdul_out[0].header.update('W_Y_END', hdul[0].header['W_Y_END'])
        hdul_out[0].header.update('NRAMPS', hdul[0].header['NRAMPS'])
        hdul_out[0].header.update('NGROUPS', hdul[0].header['NGROUPS'])
        hdul_out[0].header.update('NREADS', hdul[0].header['NREADS'])
        hdul_out[0].header.update('MJD', hdul[0].header['MJD'])
        w_f = hdul[0].header['W_Y_BEG']
        w_l = hdul[0].header['W_Y_END']
        for i in arange(1,len(hdul)):
            im = hdul[i].data - dark[w_f:w_l+1,:]
            hdul[i].data = im*coef[0] + im**2*coef[1] + \
                im**3*coef[2] + im**4*coef[3]
            hdul_out.append(fits.ImageHDU(hdul[i].data, header=hdul[i].header))
        hdul_out = fits.HDUList(hdul_out)
        hdul_out.writeto('linearity/' + f, clobber=True)    

#------------------------------------------------------------------------------#
# Frame Subtraction functions                                                  #
#------------------------------------------------------------------------------#

def frame_sub(obj):
    
    print sep
    print 'Starting single frame subtraction...'
    for f in obj.files:
        if obj.frame_mod == '2-1':
            frame_2_1(obj, f)
        elif obj.frame_mod == '3-1':
            frame_3_1(obj, f)

def frame_2_1(obj, f):
    path = obj.path
    f_out = obj.last_dir + 'mod%s_%s_%02d.fits'
    hdul = fits.open(path + f)
    print 'Processing : ', f
    ma = re.match(patt, f)
    rmps = hdul[0].header['NRAMPS']
    rds = hdul[0].header['NGROUPS']*hdul[0].header['NREADS']
    for i in range(rmps):
        hdul_out = [fits.PrimaryHDU()]
        hdul_out[0].header.update('W_Y_BEG', hdul[0].header['W_Y_BEG'], 
                              'First Y index')
        hdul_out[0].header.update('W_Y_END', hdul[0].header['W_Y_END'], 
                              'Last Y index')
        for j in range(rds-1):
            im = hdul[i*rds+j+2].data - hdul[i*rds+j+1].data
            h1 = hdul[i*rds+j+1].header
            h2 = hdul[i*rds+j+2].header
            hdu_out = fits.ImageHDU(im)
            hdu_out.header.update('EXPTIME', 
                                  around(h2['ELAPSEDT']-h1['ELAPSEDT'], 6),
                                  'Frame Exposure Time in seconds')
            hdu_out.header.update('MJD',
                                  around(hdul[0].header['MJD']+
                                         h1['ELAPSEDT']/t_c, 6),
                                  'Start of integration in MJD')
            hdul_out.append(hdu_out)

        hdul_out = fits.HDUList(hdul_out)
        write_data(hdul_out, f_out, obj, ma, i)

def frame_3_1(obj, f):
    path = obj.path
    f_out = obj.last_dir + 'mod%s_%s_%02d.fits'
    hdul = fits.open(path + f)
    print 'Processing : ', f
    ma = re.match(patt, f)
    rmps = hdul[0].header['NRAMPS']
    rds = hdul[0].header['NGROUPS']*hdul[0].header['NREADS']
    for i in range(rmps):
        hdul_out = [fits.PrimaryHDU()]
        hdul_out[0].header.update('W_Y_BEG', hdul[0].header['W_Y_BEG'], 
                              'First Y index')
        hdul_out[0].header.update('W_Y_END', hdul[0].header['W_Y_END'], 
                              'Last Y index')
        for j in range(rds-2):
            im = hdul[i*rds+j+3].data - hdul[i*rds+j+1].data
            h1 = hdul[i*rds+j+1].header
            h2 = hdul[i*rds+j+3].header
            hdu_out = fits.ImageHDU(im)
            hdu_out.header.update('EXPTIME', 
                                  around(h2['ELAPSEDT']-h1['ELAPSEDT'], 6),
                                  'Frame Exposure Time in seconds')
            hdu_out.header.update('MJD',
                                  around(hdul[0].header['MJD']+
                                         h1['ELAPSEDT']/t_c, 6),
                                  'Start of integration in MJD')
            hdul_out.append(hdu_out)

        hdul_out = fits.HDUList(hdul_out)
        write_data(hdul_out, f_out, obj, ma, i)
 
#------------------------------------------------------------------------------#
# Flat divide function                                                         #
#------------------------------------------------------------------------------#

def flat_divide(obj):
 
    flat = fits.open(obj.mast_flat)[0].data
    
    print sep
    print 'Start processing flat divided images...'
    f_in = obj.last_dir + 'mod' + obj.frame_mod + '_%s_*.fits'  
    f_out = 'flat_divide/mod%s_%s_%02d.fits'

    for f in obj.files:
        ma = re.match(patt, f)
        print 'Processing : ', f
        fs = os.popen('ls ' + f_in % ma.group(1)).read().split()
        for i in range(len(fs)):
            hdul = fits.open(fs[i])
            w_f = hdul[0].header['W_Y_BEG'] 
            w_l = hdul[0].header['W_Y_END'] 
            for j in range(1, len(hdul)):
                hdul[j].data[:,5*64:6*64] = nan
                hdul[j].data[:,30*64:] = nan
                hdul[j].data /= flat[w_f:w_l+1, :]
                if obj.flat_line:
                    flat_l = hdul[j].data
                    flat_l = median(flat_l, 0)
                    flat_l += roll(flat_l, 1) + roll(flat_l, -1) + \
                        roll(flat_l, 2) + roll(flat_l, -2)
                    flat_l = outer(zeros(w_l-w_f+1)+1, flat_l)
                    flat_l /= median(flat_l)
                    hdul[j].data /= flat_l
                hdul[j].data -= median(hdul[j].data)

            write_data(hdul, f_out, obj, ma, i)

#------------------------------------------------------------------------------#
# Pattern correction function                                                  #
#------------------------------------------------------------------------------#

def patt_corr(obj):
    
    print sep
    print 'Start pickup noise correction...'
    f_in = obj.last_dir + 'mod' + obj.frame_mod + '_%s_*.fits'  
    f_out = 'pattern_corr/mod%s_%s_%02d.fits'

    arr = arange(16)
    odd  = delete(arr, obj.odd_exc)
    evn = delete(arr, obj.evn_exc)

    for f in obj.files:
        ma = re.match(patt, f)
        print 'Processing : ', f
        fs = os.popen('ls ' + f_in % ma.group(1)).read().split()
        for i in range(len(fs)):
            hdul = fits.open(fs[i])
            w_f = hdul[0].header['W_Y_BEG'] 
            w_l = hdul[0].header['W_Y_END'] 
            y_r = w_l - w_f + 1
            for j in range(1, len(hdul)):
                im = hdul[j].data
                odd_chan = zeros([y_r, 64, odd.shape[0]])
                evn_chan = zeros([y_r, 64, evn.shape[0]])
                for k in range(odd.shape[0]):
                    odd_chan[:,:,k] = im[:,64*(2*odd[k]+1):64*(2*odd[k]+2)]
                    odd_chan[:,:,k] -= median(odd_chan[:,:,k])
                for k in range(evn.shape[0]):
                    evn_chan[:,:,k] = im[:,64*(2*evn[k]):64*(2*evn[k]+1)]
                    evn_chan[:,:,k] -= median(evn_chan[:,:,k])
                odd_corr = median(odd_chan, axis=2)
                evn_corr = median(evn_chan, axis=2)
                odd_corr -= median(odd_corr)
                evn_corr -= median(evn_corr)

                for k in arr:
                    im[:,64*(2*k+1):64*(2*k+2)] -= odd_corr
                    im[:,64*(2*k):64*(2*k+1)] -= evn_corr
                    #im[:,64*(2*k+1):64*(2*k+2)] -= \
                    #    median(im[:,64*(2*k+1):64*(2*k+2)])
                    #im[:,64*(2*k):64*(2*k+1)] -= \
                    #    median(im[:,64*(2*k):64*(2*k+1)])
                hdul[j].data = im - median(im)
               
            write_data(hdul, f_out, obj, ma, i)

#-----------------------------------------------------------------------------#
# Bad Pixels (that known a priori)                                            #
#-----------------------------------------------------------------------------#

def bad_pix_int(obj):

    mask = fits.open(obj.mask_file)[0].data

    print sep
    print 'Start processing bad pixel interpolation...'
    f_in = obj.last_dir + 'mod' + obj.frame_mod + '_%s_*.fits'  
    f_out = 'bad_pix_int/mod%s_%s_%02d.fits'

    for f in obj.files:
        ma = re.match(patt, f)
        print 'Processing : ', f
        fs = os.popen('ls ' + f_in % ma.group(1)).read().split()
        for i in range(len(fs)):
            hdul = fits.open(fs[i])
            w_f = hdul[0].header['W_Y_BEG'] 
            w_l = hdul[0].header['W_Y_END'] 
            mask_in = mask[w_f:w_l+1,:]
            ndx = where(mask_in == 1)
            for j in range(1, len(hdul)):
                ims = zeros([w_l-w_f+1, 2048, 8])
                ims[:,:,0] = roll(hdul[j].data, 1, axis=0) 
                ims[:,:,1] = roll(hdul[j].data, -1, axis=0) 
                ims[:,:,2] = roll(hdul[j].data, 1, axis=1) 
                ims[:,:,3] = roll(hdul[j].data, -1, axis=1) 
                ims[:,:,4] = roll(ims[:,:,0], 1, axis=1) 
                ims[:,:,5] = roll(ims[:,:,0], -1, axis=1) 
                ims[:,:,6] = roll(ims[:,:,1], 1, axis=1) 
                ims[:,:,7] = roll(ims[:,:,1], -1, axis=1) 
                ims = sort(ims)
                im3 = ims[:,:,3]
                im4 = ims[:,:,4]
                hdul[j].data[ndx] = (im3[ndx]+im4[ndx])*0.5


            write_data(hdul, f_out, obj, ma, i)

#----------------------------------------------------------------------------#
# Cosmic Ray Removal                                                         #
#----------------------------------------------------------------------------#

def cosmic_ray_rem(obj):

    print sep
    print 'Start processing cosmic ray removal...'
    f_in = obj.last_dir + 'mod' + obj.frame_mod + '_%s_*.fits'  
    f_out = 'cosmic_ray_rem/mod%s_%s_%02d.fits'

    for f in obj.files:
        ma = re.match(patt, f)
        print 'Processing : ', f
        fs = os.popen('ls ' + f_in % ma.group(1)).read().split()
        for i in range(len(fs)):
            hdul = fits.open(fs[i])
            w_f = hdul[0].header['W_Y_BEG'] 
            w_l = hdul[0].header['W_Y_END'] 
            mask = zeros([w_l-w_f+1, 2048])
            for j in range(1, len(hdul)):
                ims = zeros([w_l-w_f+1, 2048, 8])
                im = hdul[j].data
                ims[:,:,0] = roll(hdul[j].data, 1, axis=0) 
                ims[:,:,1] = roll(hdul[j].data, -1, axis=0) 
                ims[:,:,2] = roll(hdul[j].data, 1, axis=1) 
                ims[:,:,3] = roll(hdul[j].data, -1, axis=1) 
                ims[:,:,4] = roll(ims[:,:,0], 1, axis=1) 
                ims[:,:,5] = roll(ims[:,:,0], -1, axis=1) 
                ims[:,:,6] = roll(ims[:,:,1], 1, axis=1) 
                ims[:,:,7] = roll(ims[:,:,1], -1, axis=1) 
                ims = sort(ims)
                im3 = ims[:,:,3]
                im4 = ims[:,:,4]
                im_med = median(ims[:,:,1:7], 2) 
                im_var = var(ims[:,:,1:7], 2)
                c_m = where(abs(im - im_med) > obj.thresh*im_var)
                mask[c_m] = 1
                ndx = where(mask == 1)
                hdul[j].data[ndx] = (im3[ndx]+im4[ndx])*0.5

            write_data(hdul, f_out, obj, ma, i)

#----------------------------------------------------------------------------#
# Make a final copy                                                          #
#----------------------------------------------------------------------------#

def final_copy(obj):

    print sep
    print 'Make final copies...'
    f_in = obj.last_dir + 'mod' + obj.frame_mod + '_%s_*.fits'  
    f_out = 'reduced_data/mod%s_%s_%02d.fits'

    for f in obj.files:
        ma = re.match(patt, f)
        print 'Processing : ', f
        fs = os.popen('ls ' + f_in % ma.group(1)).read().split()
        for i in range(len(fs)):
            hdul = fits.open(fs[i])
            w_f = hdul[0].header['W_Y_BEG'] 
            w_l = hdul[0].header['W_Y_END'] 
            im = zeros([w_l-w_f+1, 2048, len(hdul)])
            for j in range(1, len(hdul)):
                im[:,:,j-1] = hdul[j].data

            if obj.comb_ramp == 1:
                hdul.append(fits.ImageHDU(median(im[:,:,:], axis=2)))
            write_data(hdul, f_out, obj, ma, i)

#-----------------------------------------------------------------------------#
# Some useful functions                                                       #
#-----------------------------------------------------------------------------#

def write_data(hdul, f_out, obj, ma, i):

    hdul.writeto(f_out % (obj.frame_mod, ma.group(1), i), clobber=True)
    print 'Ramp %02d written to %s' % (i, f_out % \
                (obj.frame_mod, ma.group(1), i))
    

#---------------------------------------------------------------------END-----#
