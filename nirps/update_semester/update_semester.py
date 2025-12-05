#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2025-04-29 at 08:48

@author: cook
"""
import os


# =============================================================================
# Define variables
# =============================================================================
INPATHS = dict()
INPATHS['user_credentials.json'] = '/scratch2/spirou/drs-bin/apero-drs-spirou-07XXX/apero/tools/resources/ari/ari-home/user_credentials.json'
INPATHS['admin.html'] = '/home/cook/ari/admin/admin.html'

# Paths to update
PATHS = dict()

PATHS['cook@jupiter'] = dict()
PATHS['nirps-client@rali'] = dict()
PATHS['spirou-client@rali'] = dict()
PATHS['cook@venus'] = dict()

# paths to sync
PATHS['cook@jupiter']['user_credentials.json'] = ['/scratch2/spirou/drs-bin/apero-drs-spirou-08XXX/apero-drs/apero/tools/resources/ari/ari-home/user_credentials.json']
PATHS['nirps-client@rali']['user_credentials.json'] = ['/cosmos99/nirps/apero-data/nirps_he_online/other/ari-home/user_credentials.json',
                                                       '/cosmos99/nirps/apero-data/nirps_ha_online/other/ari-home/user_credentials.json']
PATHS['spirou-client@rali']['user_credentials.json'] = ['/cosmos99/spirou/apero-data/spirou_offline/other/ari-home/user_credentials.json']
PATHS['cook@venus']['user_credentials.json'] = ['/export/www/home/cook/www/apero-drs/ari/home/user_credentials.json']
PATHS['cook@venus']['admin.html'] = ['/export/www/home/cook/www/apero-drs/ari/home/admin.html']

# how to sync files
SYNC = dict()
SYNC['cook@jupiter'] = 'cp {INPATH} {OUTPATH}'
SYNC['nirps-client@rali'] = 'rsync -avu {INPATH} nirps-client@rali:{OUTPATH}'
SYNC['spirou-client@rali'] = 'rsync -avu {INPATH} spirou-client@rali:{OUTPATH}'
SYNC['cook@venus'] = 'rsync -avu {INPATH} cook@venus:{OUTPATH}'


TEST = False
# -----------------------------------------------------------------------------

# =============================================================================
# Define functions
# =============================================================================
def replace_text(filename, old_text, new_text, test=False) -> bool:

    with open(filename, 'r') as file:
        content = file.read()

    if old_text not in content:
        return False

    updated_content = content.replace(old_text, new_text)

    if not test:
        with open(filename, 'w') as file:
            file.write(updated_content)

    return True

# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":
    # ----------------------------------------------------------------------
    # ask for current password
    password1 = input('Current password:\t')
    # ask for current hash
    hash1 = input('Current hash (https://asecuritysite.com/javascript/js10 HMAC Poly1305):\t')
    # ask for new password
    password2 = input('New password:\t')
    # ask for new hash
    hash2 = input('New hash (https://asecuritysite.com/javascript/js10 HMAC Poly1305):\t')

    # update local files
    for infile in INPATHS:

        print('\n\n' + '='*50)
        print(f'Updating {infile}')
        print('='*50)

        replaced1 = replace_text(INPATHS[infile], password1, password2, test=TEST)
        replaced2 = replace_text(INPATHS[infile], hash1, hash2, test=TEST)

        if replaced1:
            print(f'\nReplaced {password1}-->{password2}')
        if replaced2:
            print(f'\nReplaced {hash1}-->{hash2}')

    # update remote files
    for server in PATHS:
        print('\n\n' + '='*50)
        print(f'Updating {server}')
        print('='*50)

        for basename in PATHS[server]:

            for path in PATHS[server][basename]:

                infile = INPATHS[basename]

                command = SYNC[server].format(INPATH=infile,
                                              OUTPATH=path)

                print(f'\nRunning {command}')
                if not TEST:
                    os.system(command)

    print('\n\n\n')
    print('Don\'t forget to push the changes to github for v0.7 and v0.8')




# =============================================================================
# End of code
# =============================================================================
