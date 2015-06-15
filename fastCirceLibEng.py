#------------------------------------------------------------------------------#
# Yigit Dallilar, UF, 05/30/2015                                               #
# Contact: ydallilar@ufl.edu                                                   #
#------------------------------------------------------------------------------#
# CIRCE fast photometry pipeline tools                                         #
# -fastCirceEngine.py   version=1.0.1                                          #
# This files is the engine of the pipeline controls program flow               #
#------------------------------------------------------------------------------#
 

import os
from fastCirceLib import *
from numpy import array

allowed = ['PATH', 'FILE', 'FRAME_MOD', 'COMB_RAMP', 'BAD_PIX_MASK',
           'MASK_FILE', 'COSMIC_RAY', 'DETECT_THRESH', 'FLAT_DIVIDE', 
           'FLAT_IMAGE','FLAT_LINE', 'PATTERN_CORR', 'EVEN_EXC', 'ODD_EXC',
           'LINEARITY', 'MASTER_DARK']
sep = '-----------------------------------------------------------------------'



#-----------------------------------------------------------------------------#
# fast Circe pipeline engine class                                            #
#-----------------------------------------------------------------------------#
class fastCirceLibEng:


    def __init__(self, conf_file):

        self.conf_file = conf_file
        self.__read_config()

# read configurations

    def __read_config(self):

        print 'Reading configuration file: ', self.conf_file, ' ...'
        f = open(self.conf_file, 'r')
        confs = f.readlines()
        conf_key = []
        conf_in = []
        for i in range(len(confs)):
            if confs[i][0] != '#':
                spl = confs[i].split()
                if len(spl) == 2:
                    if if_conf_allowed(spl[0]):
                        conf_key.append(spl[0])
                        conf_in.append(spl[1])
                    else:
                        print 'WARNING: Unrecognized keyword ', \
                            spl[0], 'skipping...'
                        
        self.conf_dict = dict(zip(conf_key, conf_in))

# Initialize configurations

    def initialize(self):

        print sep
        print 'Configuration parameters:\n'

        self.last_dir = 'single_frame/'

        self.path = try_read(self,'PATH')
        error_quit('PATH', self.path)
        config_print(self, 'PATH')

        self.files = try_read(self,'FILE')
        error_quit('FILE', self.path[0])
        config_print(self, 'FILE')
        if self.files[0] == '@':
            in_fs = open(self.files[1:], 'r').read().split('\n')
            self.files = []
            for in_f in in_fs:
                if in_f.strip() != '':
                    self.files.append(in_f.strip())            
        else:
            self.files = [self.files]

        self.frame_mod = try_read(self,'FRAME_MOD')
        if self.frame_mod == -1:
            print 'No FRAME_MOD defined. Using default 2-1.'
            self.frame_mod = '2-1'
        else:
            config_print(self, 'FRAME_MOD')
            if (self.frame_mod != '2-1') & (self.frame_mod != '3-1'):
                print 'Unrecognized configuration for FRAME_MOD using 2-1 instead.'
                self.frame_mod = '2-1'

        self.comb_ramp = try_read(self,'COMB_RAMP')
        if self.comb_ramp == -1:
            print 'No COMB_RAMP defined. Enabling by default...'
            self.comb_ramp = '1'
        else:
            config_print(self, 'COMB_RAMP')
            if (self.comb_ramp != '0') & (self.comb_ramp != '1'):
                print 'Unrecognized configuration for COMB_RAMP using 1 instead.'
                self.comb_ramp = '1'
        self.comb_ramp = true_false(self.comb_ramp)
    
        self.pix_mask = try_read(self, 'BAD_PIX_MASK')
        if self.pix_mask  == -1:
            print 'No BAD_PIX_MASK defined. Disabling by default...!'
            self.pix_mask == '0'
        else:
            config_print(self, 'BAD_PIX_MASK')
            if (self.pix_mask != '0') & (self.pix_mask != '1'):
                print 'Unrecognized configuration for BAD_PIX_MASK using 0 instead.'
                self.pix_mask = '0'
            if self.pix_mask == '1':
                self.mask_file = try_read(self, 'MASK_FILE')
                if self.mask_file == -1:
                    print 'Mo MASK_FILE defined. Disabling bad pixel masking...!'
                    self.pix_mask = '0'
                else:
                    config_print(self, 'MASK_FILE')

        self.pix_mask = true_false(self.pix_mask)

        self.cosmic_ray = try_read(self, 'COSMIC_RAY')
        if self.cosmic_ray == -1:
            print 'No COSMIC_RAY defined. USing 0 instead.'
            self.cosmic_ray = 0
        else:
            config_print(self, 'COSMIC_RAY')
            if (self.cosmic_ray != '0') & (self.cosmic_ray != '1'):
                print 'Unrecognized configuration for COSMIC_RAY using 0 instead.'
                self.cosmic_ray = '0'
            if self.cosmic_ray == '1':
                self.thresh = try_read(self, 'DETECT_THRESH')
                if self.thresh == -1:
                    print 'No DETECT_THRESH defined using 5*sigma instead.'
                    self.thresh = '5'
                else:
                    config_print(self, 'DETECT_THRESH', spec=1)
                    self.thresh = float(self.thresh)
        self.cosmic_ray = true_false(self.cosmic_ray)

        self.flat_div = try_read(self, 'FLAT_DIVIDE')
        if self.flat_div == -1:
            print 'No FLAT_DIVIDE defined using 0 instead.'
            self.flat_div = '0'
        else:
            config_print(self, 'FLAT_DIVIDE')
            if (self.flat_div != '0') & (self.flat_div != '1'):
                print 'Unrecognized FLAT_DIVIDE option using 0 instead.'
                self.flat_div = '0'
            if self.flat_div == '1':
                self.mast_flat = try_read(self, 'FLAT_IMAGE')
                if self.mast_flat == -1:
                    print 'No FLAT_IMAGE defined. Disabling feature...!'
                    self.flat_div = '0'
                else:
                    config_print(self, 'FLAT_IMAGE')
        self.flat_div = true_false(self.flat_div)

        self.flat_line = try_read(self, 'FLAT_LINE')
        if self.flat_line == -1:
            print 'No FLAT_LINE defined using 0 instead'
            self.flat_line = '0'
        else:
            config_print(self, 'FLAT_LINE')
            if (self.flat_line != '1') & (self.flat_line == '0'):
                print 'Unrecognized FLAT_LINE option. Disabling feature...!'
                self.flat_line == '0'
        self.flat_line = true_false(self.flat_line)
            
        self.patt_corr_opt = try_read(self, 'PATTERN_CORR')
        if self.patt_corr_opt == -1:
            print 'No PATTERN_CORR defined using 0 instead.'
            self.patt_corr_opt == '0'
        else:
            config_print(self, 'PATTERN_CORR')
            if (self.patt_corr_opt != '1') & (self.patt_corr_opt == '0'):
                print 'Unrecognized PATTERN_CORR option. Disabling feature...!'
                self.patt_corr_opt == '0'
        self.patt_corr_opt = true_false(self.patt_corr_opt)

        self.linearity = try_read(self, 'LINEARITY')
        if self.linearity == -1:
            print 'No LINEARITY defined using 0 instead.'
            self.linearity == '0'
        else:
            config_print(self, 'LINEARITY')
            if (self.linearity != '1') & (self.linearity == '0'):
                print 'Unrecognized LINEARITY option. Disabling feature...!'
                self.linearity == '0'
            if self.linearity == '1':
                self.mst_dark = try_read(self, 'MASTER_DARK')
                if self.mst_dark == -1:
                    print 'No MASTER_DARK defined. Disabling feature...!'
                    self.mst_dark = '0'
                else:
                    config_print(self, 'MASTER_DARK')
        self.linearity = true_false(self.linearity)


        if self.patt_corr_opt == 1:
            exc = try_read(self, 'EVEN_EXC')
            if exc == -1:
                print '%-20s NO' % 'EVEN_EXC'
                self.evn_exc = array([])
            else:
                exc = array(exc.split(','))
                exc = exc.astype(int)
                out_str = '%-20s [' % 'EVEN_EXC'
                for i in exc: out_str += '%2d,' % i
                out_str = list(out_str)
                out_str[-1] = ']'
                out_str = ''.join(out_str)
                print out_str
                self.evn_exc = exc
            exc = try_read(self, 'ODD_EXC')
            if exc == -1:
                print '%-20s NO' % 'ODD_EXC'
                self.evn_exc = array([])
            else:
                exc = array(exc.split(','))
                exc = exc.astype(int)
                out_str = '%-20s [' % 'ODD_EXC'
                for i in exc: out_str += '%2d,' % i
                out_str = list(out_str)
                out_str[-1] = ']'
                out_str = ''.join(out_str)
                print out_str
                self.odd_exc = exc
                
# Talks with fastCirceLib.py

    def do_linearity(self):
        dir_to = 'linearity/'
        if not os.path.isdir(dir_to):
            os.mkdir(dir_to)
        do_linearity(self)
        self.path = 'linearity/'

    def frame_sub(self):
        dir_to = 'single_frame/'
        if not os.path.isdir(dir_to):
            os.mkdir(dir_to)
        frame_sub(self)
                    
    def flat_divide(self):
        dir_to = 'flat_divide/'
        if not os.path.isdir(dir_to):
            os.mkdir(dir_to)
        flat_divide(self)
        self.last_dir = dir_to
    
    def bad_pix_int(self):
        dir_to = 'bad_pix_int/'
        if not os.path.isdir(dir_to):
            os.mkdir(dir_to)
        bad_pix_int(self)
        self.last_dir = dir_to

    def patt_corr(self):
        dir_to = 'pattern_corr/'
        if not os.path.isdir(dir_to):
            os.mkdir(dir_to)
        patt_corr(self)
        self.last_dir = dir_to

    def cosmic_ray_rem(self):
        dir_to = 'cosmic_ray_rem/'
        if not os.path.isdir(dir_to):
            os.mkdir(dir_to)
        cosmic_ray_rem(self)
        self.last_dir = dir_to

    def final_copy(self):
        dir_to = 'reduced_data/'
        if not os.path.isdir(dir_to):
            os.mkdir(dir_to)
        final_copy(self)

# Runs all tasks together

    def runAll(self):
        if self.linearity == 1:
            self.do_linearity()
        self.frame_sub()
        if self.flat_div == 1:
            self.flat_divide()
        if self.pix_mask == 1:
            self.bad_pix_int()
        if self.patt_corr_opt == 1:
            self.patt_corr()
        if self.cosmic_ray == 1:
            self.cosmic_ray_rem()
        self.final_copy()

#-----------------------------------------------------------------------------#
# Some useful functions                                                       #
#-----------------------------------------------------------------------------#

# Check if configuration name is allowed
def if_conf_allowed(conf_key):
    check = 0  #False
    for i in range(len(allowed)):
        if conf_key == allowed[i]:
            check =  1
            break
    return check

# Can read from dict or not            
def try_read(obj, conf_key):
    try:
        info = obj.conf_dict[conf_key]
    except:
        info = -1
    return info

# Quit if any required configuration does not exist
def error_quit(conf_key, inp):
    if inp == -1:
        print 'Required paramter is not defined: "' + conf_key + '"'
        print 'Mission abort...'
        exit()

# Print option to command line        
def config_print(obj, key, spec=None):
    if spec == None:
        print "%-20s %-80s" % (key, obj.conf_dict[key])
    if spec == 1:
        print "%-20s %-80s" % (key, obj.conf_dict[key]+'*sigma')
        

# Convert to int 0 or 1   
def true_false(conf):
    if conf == '1':
        return 1
    else : 
        return 0

#-----------------------------------------------------------------------END---#
