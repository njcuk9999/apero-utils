import os
import glob
from astropy.io import fits
import numpy as np

def sync_files():
    cmd = """
    rsync -av -L --exclude={'*PCA*','*TELLU*'} -e "ssh  -oPort=5822" spirou@maestria:/space/spirou/LBL-PCA/NIRPS_HE/science  /Volumes/courlan/lbl_NIRPS_HE/
    """
    cmd = cmd.replace('\n','')

    os.system(cmd)

def create_symlinks():
    path = '/cosmos99/nirps/apero-data/nirps_he_online/out/*/*t.fits'
    path_links = '/space/spirou/LBL-PCA/NIRPS_HE/science/{}/'
    files = glob.glob(path)
    files = np.array(files)
    files = files[np.argsort(files)]

    bad_objs = ['SKY','MOON']

    for file in files:
        h = fits.getheader(file)
        if h['DRSOBJN'] in bad_objs:
            continue

        if not os.path.exists(path_links.format(h['DRSOBJN'])):
            print('creating directory {}'.format(path_links.format(h['DRSOBJN'])))
            os.mkdir(path_links.format(h['DRSOBJN']))
        cmd = 'ln -s {} {}'.format(file, path_links.format(h['DRSOBJN']))
        print(cmd)
        os.system(cmd)

    cmd = 'ln -s '
    templates = glob.glob('/cosmos99/nirps/apero-data/nirps_he_online/red/other/Template*sc1d_v*')
    template_path = '/space/spirou/LBL-PCA/NIRPS_HE/templates/'
    for i in range(len(templates)):
        template0 = templates[i]

        cmd= f'cp {template0} {template_path}'
        print(cmd)
        os.system(cmd)

def push_to_maestria():
    cmd = 'scp -r /Users/eartigau/pycodes/pixeldecorr/*.py ' \
          'spirou@maestria:/space/spirou/LBL-PCA/wraps/'
    print(cmd)
    os.system(cmd)