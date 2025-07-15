import slinky_pca
import pixpca
import wrap_lbl_slinky
import glob
import os

yaml_name = '/space/spirou/LBL-PCA/wraps/yamls/params_tellu05.yaml'

params =  pixpca.get_params(yaml_name)

slinky_pca.wrap('SPIROU')

pixpca.wrap(yaml_name)

wrap_lbl_slinky.wrapper_slinky(yaml_name)

