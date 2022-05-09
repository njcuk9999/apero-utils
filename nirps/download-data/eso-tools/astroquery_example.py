"""
Example script sent to Lison by Nicolas Buchschacher on May 1 2022.
"""
import yaml
import os
from astroquery.eso import Eso

AUTH_FILE = "auth.yaml"
with open(AUTH_FILE, 'r') as authfile:
    auth_info = yaml.safe_load(authfile)

user = auth_info["user"]
password = auth_info["password"]
password = auth_info["password"]
path_to_files = "nirps-astroquery"

eso = Eso()
eso.login(user, password)  # Ici il faut entrer le login / pass pour l’archive
eso.QUERY_INSTRUMENT_URL = "http://archive.eso.org/wdb/wdb/cas"  # Par défaut l’URL est celle de l’archive standard. Il faut donc la modifier pour NIRPS

criteria = {
    "instrument": "NIRPS",
    "night": "2022-05-04",
}  # Ici on peut mettre des critères tels que: “instrument, program, night,…"
table = eso.query_main(
    **criteria, cache=False
)  # Cette ligne fait une query sur l’archive pour lister ce qui existe
eso.retrieve_data(
    table["Dataset ID"],
    destination=path_to_files,
    unzip=False,
)  # Cette ligne lance le download et stocke les fichiers dans le répertoire mentionné plus haut.
