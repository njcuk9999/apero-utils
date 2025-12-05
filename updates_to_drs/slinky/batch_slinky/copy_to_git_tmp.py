import glob
import os

path = "/space/spirou/LBL-PCA/wraps/slinky_pca/git_temporary/apero-utils/updates_to_drs/slinky/batch_slinky"

cmd = 'rsync -arPz --include="*/" --include="*.py" --exclude="__pycache__/" --exclude="*.pyc" --exclude="*" *.py ' + path

os.system(cmd)

cmd = 'rsync -arPz yamls ' + path
os.system(cmd)

os.system('cd ' + path + ' && git add . && git commit -m "update slinky batch script" && git push')