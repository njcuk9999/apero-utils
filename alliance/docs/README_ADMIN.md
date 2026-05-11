# This is a page to remember how to do admin stuff

## How to update update apero and apero-utils

Assuming you have installed apero (with apero_install.sh)

type:
```
goapero
cat apero_instruments.conf
```

And check all instruments and repos have their correct branches.

Then to update apero and apero-utils type:
```
./apero_bin_update.sh
```

## How to update alliance apero_bin directory and docs

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


## How to install a new version of python on alliance

### Step 1

check standard software environments [here](https://docs.alliancecan.ca/wiki/Standard_software_environments)

### Step 2

Run the following:

```bash
module load StdEnv/{year} python/{py-version} hdf5
```

where `{year}` is the StdEnv year (e.g. 2023)
where `{py-version}` is the python version (e.g. 3.12)

### Step 3

Make a new virtual python environment:

```bash
cd {path}
virtualenv  {env name}
source {path}/{env name}/bin/activate
pip install --no-index --upgrade pip
```

Add or update a apero profile in:

/project/{pnumber}/apero/apero_bin/apero_profiles.conf

For example:
```
[spirou.spirou_offline_07]
module load StdEnv/2023 python/3.12 hdf5 cmake/3.31.0
source /project/6102120/apero/spirou_bin/env/apero_drs_07/bin/activate
source /project/6102120/apero/spirou_bin/settings/spirou_offline_07/spirou_offline_07.bash.setup
```

It must start with `[{INSTRUMENT}.{PROFILE_NAME}]` and then every line below
is a command that will be run when this apero profile is activated.

### Step 4 install python modules:

With v0.7.XXX:

```bash
apero-salloc --time=1:00:00 --cpus-per-task=4 --mem=8G --account=def-rdoyon
pip install -r requirements_alliance.txt
```

Note for some reason this crashes when running on the head node - do not 
run the pip install on the head node.

With v0.8.XX:

Not tested yet.

