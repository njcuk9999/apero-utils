#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2020-05-21

@author: cook
"""
import numpy as np
from tqdm import tqdm
import glob
import os
from astropy.io import fits


# =============================================================================
# Define variables
# =============================================================================
VERSION = '0.6.131'
KIND  = 'SCIENCE'
# KIND = 'TELLURIC'

if KIND == 'TELLURIC':
    IN_WORKSPACE = '/spirou/cfht_nights/cfht_sept1/reduced/'
    OUT_WORKSPACE = '/spirou/cfht_nights/telluric_targets/'
    TAR_WORKSPACE = '/home/cook/www/download/telluric_targets/'

    TARGETS = []
    TARGETS += ['74PscB','31Cas','gamTri','HR875','HR1314','pi.02Ori','HR1832',
                'zetLep','HR2180','HR2209','24Lyn','HR3131','33Lyn','etaPyx',
                '23LMi','lLeo','phiLeo','HR4687','HR4722','zetVir','82UMa',
                'HD130917','betSer','HR6025','HD159170','gamSct','51Dra','iotCyg',
                'omiCapA','chiCap','17Peg','HR8489','59Peg']

else:
    IN_WORKSPACE = '/spirou/cfht_nights/cfht_sept1/reduced/'
    OUT_WORKSPACE = '/spirou/cfht_nights/science_targets/'
    TAR_WORKSPACE = '/home/cook/www/download/science_targets/'

    TARGETS = []
    # # priority 0
    TARGETS += ['TOI-442', ' Gl699', ' Gl876', ' Gl436', ' Gl514', ' Gl382',
                ' Gl846', ' Trappist-1', ' Gl15A', ' HD189733', ' GJ1002',
                ' GJ1214']
    # priority 1
    TARGETS += ['AUMic', ' Gl388', ' TYC 4384-1735-1', ' TOI-1278', ' TOI-1759',
                ' TOI-1452']
    # priority 2
    TARGETS += ['TOI-233', ' K2-147', ' TOI-736', ' K2-33', ' TOI-876', ' TOI-732']
    # priority 3
    TARGETS += ['epsilonEri', ' GJ9827', ' Gl725A', ' Gl412A', ' Gl687', ' Gl725B']
    # priority 4
    TARGETS += ['Gl411', ' TYC3980-1081-1', ' GJ1289', ' GJ1286', ' Gl536']
    # priority 5
    TARGETS += ['V883Ori', ' DRTau', ' XZTau', ' BBW76', ' NYOri', ' V582Aur']
    # # priority 6
    TARGETS += ['TauBoo', ' V900Mon', ' V1057Cyg', ' V1515Cyg', ' WASP-127']
    # # priority 7
    TARGETS += ['HD 185603', 'RULup', 'K2-277', 'Gl905', 'HD 133112']


# extensions to look for
EXTENSIONS = ['o_pp_e2dsff_AB.fits',
              'o_pp_e2dsff_C.fits',
              'o_pp_e2dsff_tcorr_AB.fits',
              'o_pp_s1d_v_tcorr_AB.fits',
              'o_pp_s1d_v_recon_AB.fits',
              'o_pp_s1d_v_AB.fits',
              'o_pp_e2dsff_tcorr_AB_ccf_masque_sept18_andres_trans50_AB.fits',
              'o_pp_e2dsff_C_ccf_smart_fp_mask_C.fits'
              ]
# define extension header (EXTENSION[0] is reference file)
EXTHDR = 0
# odocode suffix
ODOCODE_SUFFIX = 'o_pp'
# extra files
EXTRA_FILES = ['Template_s1d_{OBJ}_sc1d_v_file_AB.fits',
               'Template_{OBJ}_tellu_obj_AB.fits']


# =============================================================================
# Define functions
# =============================================================================
def clean_objname(objname):
    objname = objname.strip()
    objname = objname.replace(' ', '_')
    objname = objname.upper()
    return objname


# =============================================================================
# Start of code
# =============================================================================
if __name__ == "__main__":

    # clean target list
    targets = clean_objname(np.char.array(TARGETS))
    # get all input files
    refpath = os.path.join(IN_WORKSPACE, '*', '2*' + EXTENSIONS[0])
    reffiles = np.sort(glob.glob(refpath))
    # storage
    objects = dict()
    # loop around reference files
    print('Group files')
    for reffile in tqdm(reffiles):
        # get odocode
        directory = os.path.dirname(reffile)
        odocode = str(os.path.basename(reffile).split(ODOCODE_SUFFIX)[0])
        # get ref header
        header = fits.getheader(reffile, ext=EXTHDR)
        # deal with header
        if 'OBJECT' in header:
            objname = clean_objname(header['OBJECT'])
        elif 'OBJNAME' in header:
            objname = clean_objname(header['OBJNAME'])
        else:
            continue
        # skip non targets
        if objname not in targets:
            continue
        # deal with new object
        if objname not in objects:
            objects[objname] = []
        # get list of files
        for ext in EXTENSIONS:
            # get filename
            filename = odocode + ext
            # get absolute path
            abspath = os.path.join(directory, filename)
            # check if it exists
            if os.path.exists(abspath):
                objects[objname].append(abspath)
        # add extra files
        kwargs = dict(OBJ=objname)
        for extrafile in EXTRA_FILES:
            # construct file name
            extrafile = extrafile.format(**kwargs)
            # locate extra files
            elistpath = os.path.join(IN_WORKSPACE, '*', extrafile)
            efilelist = np.sort(glob.glob(elistpath))
            if len(efilelist) == 0:
                continue
            # just keep the first file
            for efile in efilelist:
                if os.path.exists(efile):
                    objects[objname].append(efile)

    # make tar dir
    tar_dir = TAR_WORKSPACE
    if not os.path.exists(tar_dir):
        os.mkdir(tar_dir)
    # change to tar dir
    cwd = os.getcwd()
    os.chdir(tar_dir)
    # loop around objects
    for objname in objects:
        # print banner
        print()
        print('=' * 80)
        print('OBJECT = {0}'.format(objname))
        print('=' * 80)
        print()
        # get unique objects
        uobjects = np.unique(objects[objname])
        # make clean object name
        cobjname = objname.replace(' ', '_')
        # make output dir
        outputdir = os.path.join(OUT_WORKSPACE, cobjname)
        # remove directory if it exists
        if os.path.exists(outputdir):
            os.system('rm -rf {0}'.format(outputdir))
        # make directory
        os.mkdir(outputdir)
        # make a sim link in this dir
        print('Making symlinks (OBJECT={0})'.format(objname))
        for filename in tqdm(uobjects):
            # make target file path
            dst = os.path.join(outputdir, os.path.basename(filename))
            # make sim link
            os.symlink(filename, dst)
        # tar dir
        print('Making tar.gz for OBJECT={0}'.format(objname))
        dst = '{0}_v{1}.tar.gz'.format(cobjname, VERSION)
        command = 'tar -hczvf {0} {1}'.format(dst, outputdir)
        os.system(command)
    # change back to current dir
    os.chdir(cwd)

# =============================================================================
# End of code
# =============================================================================
