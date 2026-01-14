# Temporary babysitter instructions


## How to install apero for a new user


This only needs to be done once per instrument:

type:
```bash
/project/612120/apero/apero_bin/apero_install.sh {instrument}
```

If you leave the instrument blank you should see a list of available instruments.

If you have installed this instrument before this script will do nothing.


## How to Activate the correct apero profile

Activate in the head node (the one yuo logged into)

Command syntax is as follows:
```
apero-activate {instrument} {apero profile}
```

This will list the instruments available
```
apero-activate 
```

This will list the apero profiles available for spirou
```
apero-activate spirou
```

Once you see the following it means you are ready to continue:

```
HH:MM:SS.SS-**|VALID| ***************************************************************************
HH:MM:SS.SS-**|VALID| Recipe apero_validate has been successfully completed	(N.NNN seconds)
HH:MM:SS.SS-**|VALID| ***************************************************************************
```


## How to run the checks

In the head node (the ony you logged into) type the following commands for:

- raw checks:
```
apero-checks
python apero_raw_check.py {check yaml} {other args}
```

- red checks:
```
apero-checks
python apero_red_check.py {check yaml} {other args}
```

Where 'other args' are usually one of the following:
```
--today
--yesterday
--obsdir={DATE}
```
And if wanting test mode using the --test={TEST NAME} argument.


If you have not activate an apero profile this will not work!


## How to run the manual trigger

In the head node (the one you logged into) type the following commands:

```
apero-trigger
python manual_trigger.py {trigger yaml} --obsdir={DATE}
```

Note that the trigger yaml should be as follows:
- SPIROU: spirou_profile_online_alliance.yaml
- NIRPS-HE: nirps_he_profile_online_alliance.yaml
- NIRPS-HA: nirps_ha_profile_online_alliance.yaml

In general it should follow {apero profile}_alliance.yaml syntax

If you have not	activate an apero profile this will not	work!
