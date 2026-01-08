#!/bin/bash

# Note this code is only to be run once per user
#    Do not add anything here to change setups
#    This just puts stuff in the ~/.bashrc
#    Please look at apero_instruement.ini to see list of setup scripts

# Get core source script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $SCRIPT_DIR/apero_core.sh

# Find the bash file to push this into
TARGET="$HOME/.bashrc"

# Set the path to profiles.ini
INSTRUMENT_FILE="$APERO_BIN_PATH/apero_instruments.ini"


# -----------------------------------------------------------------------------
# Function: Show help message
# -----------------------------------------------------------------------------
show_help() {
    echo ""
    echo "Usage: apero_install.sh <instrument_name> [--debug]"
    echo ""
    echo "This script should NOT be sourced."
    echo ""
    echo "This script sets up the specified instrument profile by adding the"
    echo "necessary source commands to your ~/.bashrc (only once per user)."
    echo ""
    echo "Available instruments:"
    if [[ -f "$INSTRUMENT_FILE" ]]; then
        grep -o '^[^=]*' "$INSTRUMENT_FILE" | sed 's/^/  - /'
    else
        echo "  (Instrument file not found: $INSTRUMENT_FILE)"
    fi
    echo ""
    echo "Example:"
    echo "  ./apero_install.sh nirps"
    echo "  ./apero_install.sh spirou --debug"
    echo ""
    echo "Notes:"
    echo "  - The --debug flag prints internal paths for troubleshooting."
    echo "  - Do not run this script multiple times; it will append only once."
    echo ""
}
# -----------------------------------------------------------------------------
# Detect if script is sourced, and exit if it is
# -----------------------------------------------------------------------------
if [ "${BASH_SOURCE[0]}" != "$0" ]; then
    echo "ERROR: This script should NOT be sourced."
    show_help
    return 1 2>/dev/null || exit 1
fi
# -----------------------------------------------------------------------------
#  Check for help flag
# -----------------------------------------------------------------------------
if [[ "$1" == "-h" || "$1" == "--help" ]]; then
    show_help
    exit 0
fi

# -----------------------------------------------------------------------------
# Check if apero_instruments.ini exists
# -----------------------------------------------------------------------------
if [ ! -f "$INSTRUMENT_FILE" ]; then
  echo "Error: $INSTRUMENT_FILE file not found."
  exit 1
fi

# -----------------------------------------------------------------------------
# Check if instrument name is provided as an argument
# -----------------------------------------------------------------------------
# Check if instrument name is provided as an argument
if [ -z "$1" ]; then
    show_help
    exit 1
fi

# -----------------------------------------------------------------------------
# print that we found it (may remove later)
# -----------------------------------------------------------------------------
if [ "$2" = "--debug" ]; then
  echo "Found: INSTRUMENT_FILE=$INSTRUMENT_FILE"
fi

# -----------------------------------------------------------------------------
# Read profiles.ini and find the path for the provided profile
# -----------------------------------------------------------------------------
# Check if profile path exists
INSTRUMENT_PATH=$(grep "^$1=" "$INSTRUMENT_FILE" | cut -d'=' -f2)
if [ -z "$INSTRUMENT_PATH" ]; then
    echo "Error: Instrument '$1' not found in $INSTRUMENT_FILE."
    show_help
    exit 1
fi

# -----------------------------------------------------------------------------
# print that we found it (may remove later)
# -----------------------------------------------------------------------------
if [ "$2" = "--debug" ]; then
  echo "Found: INSTRUMENT_PATH=$INSTRUMENT_PATH"
fi

# -----------------------------------------------------------------------------
# Determine the setup script to run based on OS
# -----------------------------------------------------------------------------
SETUP_SCRIPT="$APERO_BIN_PATH/$INSTRUMENT_PATH"

# -----------------------------------------------------------------------------
# Check if the setup script file actually exists
# -----------------------------------------------------------------------------
if [ ! -f "$SETUP_SCRIPT" ]; then
  echo ""
  echo "Error: Setup script for instrument '$1' not found at:"
  echo "  $SETUP_SCRIPT"
  echo ""
  echo "The instrument is defined in $INSTRUMENT_FILE but the script file is missing."
  echo ""
  echo "Available instruments in $INSTRUMENT_FILE:"
  grep -o '^[^=]*' "$INSTRUMENT_FILE"
  echo ""
  exit 1
fi

read -r -d '' SNIPPET <<EOF
# Source $1 profile if present
if [ -f "$SETUP_SCRIPT" ]; then
    source "$SETUP_SCRIPT"
fi
EOF

# -----------------------------------------------------------------------------
# print that we found it (may remove later)
# -----------------------------------------------------------------------------
if [ "$2" = "--debug" ]; then
  echo "SNIPPET = $SNIPPET"
fi

# -----------------------------------------------------------------------------
# Append only if not present
# -----------------------------------------------------------------------------
if ! grep -Fq "$SETUP_SCRIPT" "$TARGET"; then
    echo ""
    echo "Installing instrument profile for '$1' into $TARGET"
    echo ""
    printf "\n%s\n" "$SNIPPET" >> "$TARGET"
    source "$SETUP_SCRIPT"
    source ~/.bashrc
else
    if [ -f "$SETUP_SCRIPT" ]; then
        if [ "$2" = "--debug" ]; then
          echo "Not installing (already present in $TARGET)"
        fi
        source "$SETUP_SCRIPT"
    else
        if [ "$2" = "--debug" ]; then
          echo "Not installing (already present in $TARGET) and not sourcing"
        fi
    fi
fi




