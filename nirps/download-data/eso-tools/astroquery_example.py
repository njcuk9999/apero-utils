"""
Example script sent to Lison by Nicolas Buchschacher on May 1 2022.
"""
import os
from astroquery.eso import Eso

user = "TODO"
password = "TODO"
path_to_files = "TODO"

eso = Eso()
eso.login(user, password)  # Ici il faut entrer le login / pass pour l’archive
eso.QUERY_INSTRUMENT_URL = "http://archive.eso.org/wdb/wdb/cas"  # Par défaut l’URL est celle de l’archive standard. Il faut donc la modifier pour NIRPS
os.environ["XDG_CACHE_HOME"] = [
    path_to_files
]  # Ici tu as la possibilité de contrôler où seront stockés les fichiers téléchargés

criteria = {
    "instrument": "NIRPS",
    "night": "2022-04-26",
}  # Ici on peut mettre des critères tels que: “instrument, program, night,…"
table = eso.query_main(
    **criteria, cache=False
)  # Cette ligne fait une query sur l’archive pour lister ce qui existe
eso.retrieve_data(
    table["Dataset ID"]
)  # Cette ligne lance le download et stocke les fichiers dans le répertoire mentionné plus haut.
