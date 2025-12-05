
import glob 
import os 



fichiers_a_verifier = [
    '2412942o.fits', '2412943o.fits', '2412944o.fits', '2412945o.fits',
    '2412946o.fits', '2412947o.fits', '2412948o.fits', '2412949o.fits',
    '2815524o.fits', '2815524o.fits', '2839133o.fits', '2839134o.fits',
    '2839133o.fits', '2839134o.fits', '2815396o.fits', '2815396o.fits',
    '2815397o.fits', '2815397o.fits', '2838840o.fits', '2838840o.fits',
    '2838841o.fits', '2838841o.fits', '2838997o.fits', '2838997o.fits',
    '2838998o.fits', '2838998o.fits', '2838999o.fits', '2838999o.fits',
    '2839000o.fits', '2839000o.fits', '3029164o.fits', '3029164o.fits',
    '3029346o.fits', '3029346o.fits', '2412950o.fits', '2412950o.fits',
    '2412951o.fits', '2412951o.fits', '2412952o.fits', '2412952o.fits',
    '2412953o.fits', '2412953o.fits', '2412954o.fits', '2412954o.fits',
    '2412955o.fits', '2412955o.fits', '2412956o.fits', '2412956o.fits',
    '2412957o.fits', '2412957o.fits', '2812872o.fits', '2812872o.fits',
    '2812873o.fits', '2812873o.fits', '2812874o.fits', '2812874o.fits'
]

nights = []
full_path = []

# Dossier de base (optionnel, ici c'est le dossier courant) 
base_dir = '.' 

# Pour chaque fichier, on cherche dans les sous-dossiers 
for nom_fichier in fichiers_a_verifier: 
    # Génère le pattern de recherche 
    pattern = os.path.join(base_dir, '**', nom_fichier) 
    # Recherche récursive 
    correspondances = glob.glob(pattern, recursive=True) 
    # Filtrage : on ne garde que les fichiers valides 
    fichiers_existants = [f for f in correspondances if os.path.isfile(f)] 
    if fichiers_existants: 
        print(f"✅ Fichier '{nom_fichier}' trouvé dans : {fichiers_existants}") 

        night  = fichiers_existants[0].split('/')[1]
        nights.append(night)
        full_path.append(os.path.abspath(fichiers_existants[0]))
    else: 
        print(f"❌ Fichier '{nom_fichier}' non trouvé dans les sous-dossiers.") 
        nights.append('None')
        full_path.append('None')


path0 = '/cosmos99/spirou/apero-data/spirou_offline/tmp/{}/{}*'
is_missing = []

for ifile in range(len(fichiers_a_verifier)):
    fichier = fichiers_a_verifier[ifile]
    night = nights[ifile]
    path = path0.format(night, fichier.split('.')[0])

    # Recherche récursive 
    correspondances = glob.glob(path, recursive=True) 
    
    if len(correspondances) > 0:
        # Filtrage : on ne garde que les fichiers valides
        print(f"✅ Dossier '{path}' trouvé.")
    else:
        print(f"❌ Dossier '{path}' non trouvé.")
        is_missing.append(full_path[ifile])

# we get the headers for the missing files
from astropy.io import fits

object_missing = []
nights_missing = []
for ifile in range(len(is_missing)):
    fichier = is_missing[ifile]

    hdr = fits.getheader(fichier)
    obj = hdr['OBJECT']
    if obj not in object_missing:
        object_missing.append(obj)
    if nights[ifile] not in nights_missing:
        nights_missing.append(nights[ifile]) 

    print(obj, nights[ifile], full_path[ifile])

print('Objects missing:')
for obj in object_missing:
    print(obj)
print('Nights missing:')
for night in nights_missing:
    print(night)