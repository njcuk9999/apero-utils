# Set the project ID (it may change in future)
export APERO_PROJECT_ID="6102120"
export APERO_SERVER="alliance"

# define the apero bin path
APERO_BIN_PATH="/project/$APERO_PROJECT_ID/apero/apero_bin"

# global location aliases
alias goapero="cd /project/$APERO_PROJECT_ID/apero/"

# software aliases
alias dfits="/project/$APERO_PROJECT_ID/apero/spirou_bin/scripts/fitsio/dfits"
alias fitsort="/project/$APERO_PROJECT_ID/apero/spirou_bin/scripts/fitsio/fitsort"

# apero tools
alias apero-trigger="cd /project/6102120/apero/spirou_bin/scripts/apero-utils/nirps/manual_trigger"
alias apero-checks="cd /project/6102120/apero/spirou_bin/scripts/apero-utils/nirps/apero_check"

# activate apero profiles
#   please add apero profiles to apero_profiles.conf
alias apero-activate="source $APERO_BIN_PATH/apero_activate.sh"
