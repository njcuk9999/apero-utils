#!/bin/bash

# Get core source script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $SCRIPT_DIR/apero_core.sh

# Set instrument
APERO_INSTRUMENT=$1

# -----------------------------------------------------------------------------
# Function: Show help for APERO environment setup
# -----------------------------------------------------------------------------
show_help() {
    echo ""
    echo "APERO Environment Setup [$APERO_INSTRUMENT]"
    echo "----------------------"
    echo
    echo "Please only source this file"
    echo
    echo "Project ID:      $APERO_PROJECT_ID"
    echo "Server:          $APERO_SERVER"
    echo "Instrument:      $APERO_INSTRUMENT"
    echo "APERO Path:      $APERO_PATH"
    echo "APERO Bin Path:  $APERO_BIN_PATH"
    echo ""
    echo "Available aliases:"
    echo "  goapero        : cd to the APERO project directory"
    echo "  gobin          : cd to the $APERO_INSTRUMENT bin directory"
    echo "  godata         : cd to the $APERO_INSTRUMENT data directory"
    echo "  dfits          : run dfits for $APERO_INSTRUMENT"
    echo "  fitsort        : run fitsort for $APERO_INSTRUMENT"
    echo "  apero-trigger  : cd to manual trigger scripts for $APERO_INSTRUMENT"
    echo "  apero-checks   : cd to APERO checks for $APERO_INSTRUMENT"
    echo "  apero-activate : source the APERO profile activation script"
    echo "  apero-salloc   : run APERO salloc launcher"
    echo "  apero-find     : run APERO file finder tool"
    echo ""
    echo "Usage example:"
    echo "  source this_script.sh      # sets up environment and aliases"
    echo "  goapero                    # quickly navigate to APERO directory"
    echo "  dfits file.fits            # run dfits on a FITS file"
    echo ""
}

# -----------------------------------------------------------------------------
#  Check for help flag
# -----------------------------------------------------------------------------
if [[ "$1" == "-h" || "$1" == "--help" ]]; then
    APERO_INSTRUMENT="NOT SET"
    show_help
    return 0
fi

# -----------------------------------------------------------------------------
# Functions to change directory to the instrument bin and data directory
# -----------------------------------------------------------------------------
gofunc() {
    # $1 = type (bin, data, etc.)
    # $2 = instrument

    local type="$1"
    local inst="$2"

    # deal with no argument
    if [[ -z "$inst" ]]; then
        echo "Usage: go${type} <instrument>"
        echo ""
        echo "Change to the <instrument> ${type} directory."
        return 1
    fi

    # help
    if [[ "$inst" == "-h" || "$inst" == "--help" ]]; then
        echo "Usage: go${type} <instrument>"
        echo ""
        echo "Change to the <instrument> ${type} directory."
        return 0
    fi

    # set path
    local path="$APERO_PATH/${inst}_${type}"

    # check directory
    if [[ ! -d "$path" ]]; then
        echo "*************************"
        echo "ERROR: directory does not exist:"
        echo "  $path"
        echo "*************************"
        return 1
    fi

    cd "$path" || return
}

gobin() {
    gofunc "bin" "$@"
}
godata() {
    gofunc "data" "$@"
}

# -----------------------------------------------------------------------------
#  Set global variables and aliases
# -----------------------------------------------------------------------------
# global location aliases
alias goapero="cd $APERO_PATH"

# activate apero profiles
#   please add apero profiles to apero_profiles.conf
alias apero-activate="source $APERO_BIN_PATH/apero_activate.sh"

# apero salloc launcher
alias apero-salloc="$APERO_BIN_PATH/apero_salloc.sh"

# apero find launcher
alias apero-find="$APERO_BIN_PATH/apero_find.sh"

# -----------------------------------------------------------------------------
#  Detect if not sourced - show help and exit
# -----------------------------------------------------------------------------
(return 0 2>/dev/null)
if [ $? -ne 0 ]; then
    show_help
    exit 1
fi

