#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Copy a ARI list or directory of files (or paths copied via user input)

i.e.

python ari_smart_download.py Delgado_Mena_2023_07_26/files
       --dir=Delgado_Mena_2023_07_26/lists --tar=Delgado_Mena_2023_07_26.tar.gz

will produce a tar file of all files in the "lists" directory

Created on 2023-07-26

@author: cook, artigau
"""
import os
import numpy as np
import argparse
import tarfile
import glob
from astropy.io import fits

# =============================================================================
# Define variables
# =============================================================================
# define SYNC type (cp, scp, rsync)
SYNC_TYPE = 'rsync'
# this is the command to use with SYNC_TYPE
if SYNC_TYPE == 'cp':
    SYNC_CMD = 'cp {INPATH} {OUTPATH}'
elif SYNC_TYPE == 'scp':
    SYNC_CMD = 'scp -r {USER}@{HOST}:{INPATH} {OUTPATH}'
elif SYNC_TYPE == 'rsync':
    SYNC_CMD = 'rsync -avu {USER}@{HOST}:{INPATH} {OUTPATH}'
else:
    SYNC_CMD = 'cp {INPATH} {OUTPATH}'

# fill these out if you want to scp/rsync: note you must have tunnelling
# enabled to do this. Set to None to just copy to a path
HOST = 'rali'
USER = 'nirps-client'

# =============================================================================
# define functions
# =============================================================================
def download(infilename, outfilename):
    # set up dictionary
    cmd_kwargs = dict()
    cmd_kwargs['INPATH'] = infilename
    cmd_kwargs['OUTPATH'] = outfilename
    cmd_kwargs['HOST'] = HOST
    cmd_kwargs['USER'] = USER
    cmd = SYNC_CMD.format(**cmd_kwargs)
    # try to copy / download file
    try:
        print(cmd)
        if not args.test:
            os.system(cmd)
    except Exception as e:
        emsg = 'Cannot run {0} \n\t Error {1}: {2}'
        eargs = [cmd, type(e), str(e)]
        print(emsg.format(*eargs))


# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":
    # set up argparse
    parser = argparse.ArgumentParser(description='Smart download of ARI files')
    parser.add_argument('path', type=str, action='store',
                        help='The local path to save files to (if using --tar'
                             ' this should be an empty directory)')
    parser.add_argument('--dir', type=str, action='store', default=None,
                        help='A directory containing only ARI lists')
    parser.add_argument('--files', type=list, action='store', nargs='+',
                        help='A list of files containing paths from ARI, if '
                             'blank user is asked to copy and paste files')
    parser.add_argument('--tar', action='store', type=str,
                        default=None,
                        help='Tar up all products in "path" into this tarname')
    parser.add_argument('--test', action='store_true',
                        help='Test mode (no file operations)')
    parser.add_argument('--seedpath', action='store', type=str, default=None,
                        help='If using ARI remove this path from the start '
                             'of files and keep the rest of the file structure'
                             ' the same. Should be a block path '
                             '(e.g. end in tmp, red, calib, tellu, out etc)')
    parser.add_argument('--get_calibs', action='store_true', default=False,
                        help='Get matching calibration files for the ARI files.'
                             'Note requires --calibpath and --seedpath to be '
                             'set.')
    parser.add_argument('--calibpath', action='store', type=str, default=None,
                        help='The path on the server where calibrations are '
                             'stored.')
    parser.add_argument('--localcalibpath', action='store', type=str,
                        default=None,
                        help='The path on the local machine where calibrations '
                             'should be saved to.')
    parser.add_argument('--get_tellurics', action='store_true', default=False,
                        help='Get matching telluric files for the ARI files.')
    parser.add_argument('--tellupath', action='store', type=str, default=None,
                        help='The path on the server where tellurics are '
                             'stored.')
    parser.add_argument('--localtellupath', action='store', type=str,
                        default=None,
                        help='The path on the local machine where tellurics '
                             'should be saved to.')
    # get arguments
    args = parser.parse_args()
    # -------------------------------------------------------------------------
    # deal with path not existing
    if not os.path.exists(args.path):
        os.makedirs(args.path)
    # -------------------------------------------------------------------------
    # deal with being giving a directory
    if args.dir is not None:
        if os.path.exists(args.dir):
            args.files = glob.glob(os.path.join(args.dir, '*.txt'))
    # -------------------------------------------------------------------------
    # if we have files then get paths from ARI lists
    if args.files is not None:
        # empty list to store user inputs
        user_inputs = dict()
        # loop around files
        for filename in args.files:
            # only deal with text files
            if not filename.endswith('.txt'):
                continue
            # get basename
            basename = os.path.basename(filename).replace('.txt', '')
            if os.path.exists(filename):
                # read lines
                with open(filename, 'r') as pathfile:
                    lines = pathfile.readlines()
                # need to remove escape characters
                valid_files = []
                for line in lines:
                    valid_file = line.strip('\n')
                    if len(valid_file) > 0:
                        valid_files.append(valid_file)
                # push lines into user_inputs
                user_inputs[basename] = valid_files
    else:
        # Copy-paste the name of all files for the ARI interface "file_list" document.
        # Press enter twice when done and files get copied to the current directory.
        # Note that unless you have an ssh key, you will be prompted for your password
        # for each file. This is a feature, not a bug. Look at this website for more
        # information on how to set these keys: https://www.ssh.com/ssh/copy-id

        # empty list to store user inputs
        user_inputs = dict()
        user_inputs[''] = []
        user = ' '  # initialize user input
        print('Copy-paste the name of all files. Press enter twice when done.\n')
        while user != '':
            user = input()
            user_inputs[''].append(user)
    # -------------------------------------------------------------------------
    # save a list of files
    outpaths = []
    # loop around user_inputs
    for user_input in user_inputs:
        # if user inputs is blank we do not want a subdir
        if len(user_input) == 0:
            outpath = args.path
        # otherwise we put outputs in directory
        else:
            outpath = os.path.join(args.path, user_input)
        # make sure outpath exists locally
        if not os.path.exists(outpath):
            os.makedirs(outpath)
        # we now copy files
        for it, infilename in enumerate(user_inputs[user_input]):
            # skip blank file names
            if len(infilename) == 0:
                continue
            # deal with seedpath
            if args.seedpath is not None:
                # get the seedpath
                seedpath = args.seedpath
                # get the infilename
                outfilename = os.path.abspath(infilename)
                # remove the seedpath from the infilename
                outfilename = outfilename.replace(seedpath, '')
                # remove the first character (should be a /)
                if outfilename[0] == os.sep:
                    outfilename = outfilename[len(os.sep):]
                # add the outpath to the outfilename
                outfilename = os.path.join(outpath, outfilename)
                # create absolute outpath for filename
                if not os.path.exists(os.path.dirname(outfilename)):
                    os.makedirs(os.path.dirname(outfilename))
            else:
                # create absolute outpath for filename
                outfilename = os.path.join(outpath, os.path.basename(infilename))
            # -----------------------------------------------------------------
            # store outfiles
            outpaths.append(outfilename)
            # -----------------------------------------------------------------
            # skip files we already have
            if os.path.exists(outfilename):
                print('File {} already exists. '
                      'Skipping download.'.format(outfilename))
                continue
            # if we are copying locally we can check that the inpath exist
            if SYNC_TYPE == 'cp':
                if not os.path.exists(infilename):
                    print('File {} does not exist.'
                          'Skipping download'.format(infilename))
                    continue
            # -----------------------------------------------------------------
            # print progress
            print('*' * 50)
            margs = [infilename, it + 1, len(user_inputs[user_input])]
            print('Getting {0} [{1}/{2}]'.format(*margs))
            print('*' * 50)
            # download the files
            download(infilename, outfilename)

    # -------------------------------------------------------------------------
    # get calibrations if required
    if args.get_calibs:
        if args.localcalibpath is None:
            raise ValueError('Must set --localcalibpath to get calibration files')
        elif not os.path.exists(args.localcalibpath):
            os.makedirs(args.localcalibpath)
        if args.calibpath is None:
            raise ValueError('Must set --calibpath to get calibration files')

        # loop around outpaths
        for it, outfile in enumerate(outpaths):
            # print progress
            print('*' * 50)
            margs = [outfile, it + 1, len(outpaths)]
            print('Getting calibrations for {0} [{1}/{2}]'.format(*margs))
            print('*' * 50)
            # read the header
            header = fits.getheader(outfile)
            # get all the calib file keys
            calib_keys = header['CDB*']
            # loop around calib files
            for calib_key in calib_keys:
                # get the calib file
                calib_file = header[calib_key]
                # deal with no file (e.g. "No leak")
                if calib_file.upper().startswith('NO '):
                    continue
                # get the calib file name
                calib_filename = os.path.join(args.calibpath, calib_file)
                # get the local calib file name
                local_calib_filename = os.path.join(args.localcalibpath,
                                                    calib_file)
                # deal with file already on disk
                if os.path.exists(local_calib_filename):
                    print('File {} already exists. '
                          'Skipping download.'.format(local_calib_filename))
                    continue
                # download the files
                download(calib_filename, local_calib_filename)
            # delete header
            del header
    # -------------------------------------------------------------------------
    # get telluric files if required
    if args.get_tellurics:
        # print progress
        print('*' * 50)
        print('Getting telluric files')
        print('*' * 50)
        # deal with no tellupath
        if args.tellupath is None:
            raise ValueError('Must set --tellupath to get telluric files')
        # deal with no localtellupath
        if args.localtellupath is None:
            raise ValueError('Must set --localtellupath to get telluric files')
        elif not os.path.exists(args.localtellupath):
            os.makedirs(args.localtellupath)
        # download the trans models
        trans_model = 'trans_model_*.fits'
        download(os.path.join(args.tellupath, trans_model), args.localtellupath)
        # download the trans files
        tellu_files = '*pp_tellu_trans_*.fits'
        download(os.path.join(args.tellupath, tellu_files), args.localtellupath)
    # -------------------------------------------------------------------------
    # tar up all products if required
    if args.tar is not None:
        # get tar directory
        tardir = os.path.dirname(args.tar)
        if len(tardir) == 0:
            tardir = os.getcwd()
        # deal with tar directory not existing
        if not os.path.exists(tardir):
            os.makedirs(tardir)

        print('Compressing {0} to {1}'.format(args.path, args.tar))
        if not args.test:
            # add to the archive
            with tarfile.open(args.tar, 'w:gz') as tar:
                tar.add(args.path)
    # -------------------------------------------------------------------------
    print('Code finished successfully')



