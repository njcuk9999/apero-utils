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

# =============================================================================
# Define variables
# =============================================================================
# define SYNC type (cp, scp, rsync)
SYNC_TYPE = 'cp'
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
HOST = None
USER = 'nirps-client'

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
        i = 1  # initialize counter
        print('Copy-paste the name of all files. Press enter twice when done.\n')
        while user != '':
            user = input()
            user_inputs[''].append(user)
            i += 1
        # make a numpy array
        user_inputs = np.array(user_inputs)

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
        for infilename in user_inputs[user_input]:
            # create absolute outpath for filename
            outfilename = os.path.join(outpath, os.path.basename(infilename))
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
                continue
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



