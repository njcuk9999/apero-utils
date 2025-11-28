# Set the project ID (it may change in future)
export APERO_PROJECT_ID="6102120"

# define the apero bin path
APERO_BIN_PATH="/project/$APERO_PROJECT_ID/apero/apero_bin"

# global location aliases
alias goapero="cd /project/$APERO_PROJECT_ID/apero/"

# software aliases
alias dfits="/project/$APERO_PROJECT_ID/apero/spirou_bin/scripts/fitsio/dfits"
alias fitsort="/project/$APERO_PROJECT_ID/apero/spirou_bin/scripts/fitsio/fitsort"

# activate apero profiles
#   please add apero profiles to apero_profiles.conf
alias apero-activate="source $APERO_BIN_PATH/apero_activate.sh"
