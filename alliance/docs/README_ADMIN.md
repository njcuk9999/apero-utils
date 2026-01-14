# This is a page to remember how to do admin stuff





## How to update alliance apero_bin directory

Assuming you have installed apero (with apero_install.sh)

type:
```
goapero
apero_bin/apero_bin_update.sh
```



## How to add aliases to instrument scripts 


From github get apero-utils


Edit `apero-utils/alliance/apero_bin/apero_activate.sh`


Specifically look for the section # software aliases

Commit your changes.
Note the apero-utils repos must be on the branch you commited to for changes to occur.
On Alliance please check:
- For SPIROU: `/project/6102120/apero/spirou_bin/scripts/apero-utils`
- For NIRPS: `/project/6102120/apero/nirps_bin/scripts/apero-utils`

with `git branch`

On Alliance (on the head node)

Assuming you've installed apero (with apero_install.sh)

type:
```
goapero
apero_bin/apero_bin_update.sh
```

This should update the scripts for everyone.
