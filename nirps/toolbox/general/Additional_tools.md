# Additional tools

author: Neil Cook
last update: 2023-03-13

## Apero reduction interface (ARI)

URL to visit page: http://apero.exoplanets.ca/ari/

    cd /cosmos99/spirou/apero-bin/apero-utils/general/apero_reduction_interface
    python simple_ari.py {profile name}.yaml

where:

- `{profile name}` is one of the yamls in the `apero_reduction_interface` directory
  - These list the apero profiles and directories for each reduction
  - Note all previous reductions will be removed from the ARI tables use the skip options to not regenerate 
other profiles which you do not wish to update now


## Apero manual trigger

Current procedure:

- make sym links
- apero_precheck [not functional yet]
- apero_processing
- apero_get 
- ARI [not functional yet]
- LBL [not functional yet]
  - requires matching templates to objects


    cd /cosmos99/spirou/apero-bin/apero-utils/nirps/manual_trigger
    python manual_trigger.py {profile name}.yaml --obsdir={obs dirs}

where:

- `{profile name}` is one of the yamls in the `manual_trigger` directory
- `{obs dirs}` is a single observation directory or a comma no white space separated list of directories
  - i.e. 2023-03-10
  - i.e. 2023-03-10,2023-03-11
  - Note you can leave out the `--obsdir` to reduce all directories
  - Note these 

