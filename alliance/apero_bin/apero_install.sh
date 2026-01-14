#!/bin/bash

# Note this code is only to be run once per user
#    Do not add anything here to change setups
#    This just puts stuff in the ~/.bashrc
#    Please look at apero_instruments.conf to see list of setup scripts

# Get core source script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $SCRIPT_DIR/apero_core.sh

# Find the bash file to push this into
TARGET="$HOME/.bashrc"

# Set the path to instruments conf
INSTRUMENT_FILE="$APERO_BIN_PATH/apero_instruments.conf"


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
        # list section names [instrument]
        grep -E '^\[[^]]+\]' "$INSTRUMENT_FILE" | sed 's/^\[/- /; s/\]$//' | sed 's/^/  /'
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
# Check if apero_instruments.conf exists
# -----------------------------------------------------------------------------
if [ ! -f "$INSTRUMENT_FILE" ]; then
  echo "Error: $INSTRUMENT_FILE file not found."
  exit 1
fi

# -----------------------------------------------------------------------------
# Check if instrument name is provided as an argument
# -----------------------------------------------------------------------------
if [ -z "$1" ]; then
    show_help
    exit 1
fi
INSTRUMENT_NAME="$1"

# -----------------------------------------------------------------------------
# print that we found it (may remove later)
# -----------------------------------------------------------------------------
if [ "$2" = "--debug" ]; then
  echo "Found: INSTRUMENT_FILE=$INSTRUMENT_FILE"
fi

# -----------------------------------------------------------------------------
# Resolve install_script for the provided instrument from the conf file
# -----------------------------------------------------------------------------
# parse the file: track when inside [instrument], then read install_script
INSTALL_SCRIPT_REL=""
CURRENT_SECTION=""
while IFS= read -r line || [ -n "$line" ]; do
    # trim whitespace
    line="$(echo "$line" | sed 's/^\s\+//;s/\s\+$//')"
    # skip comments/empty
    [[ -z "$line" || "$line" =~ ^# ]] && continue

    if [[ "$line" =~ ^\[(.+)\]$ ]]; then
        CURRENT_SECTION="${BASH_REMATCH[1]}"
        continue
    fi
    if [[ "$CURRENT_SECTION" == "$INSTRUMENT_NAME" ]]; then
        if [[ "$line" =~ ^install_script=(.+)$ ]]; then
            INSTALL_SCRIPT_REL="${BASH_REMATCH[1]}"
            break
        fi
    fi

done < "$INSTRUMENT_FILE"

if [[ -z "$INSTALL_SCRIPT_REL" ]]; then
    echo "Error: Instrument '$INSTRUMENT_NAME' not found or missing install_script in $INSTRUMENT_FILE."
    echo "Available instruments:"
    grep -E '^\[[^]]+\]' "$INSTRUMENT_FILE" | sed 's/^\[//; s/\]$//' | sed 's/^/  - /'
    exit 1
fi

# -----------------------------------------------------------------------------
# print that we found it (may remove later)
# -----------------------------------------------------------------------------
if [ "$2" = "--debug" ]; then
  echo "Found: INSTALL_SCRIPT_REL=$INSTALL_SCRIPT_REL"
fi

# -----------------------------------------------------------------------------
# Determine the setup script to run based on conf
# -----------------------------------------------------------------------------
SETUP_SCRIPT="$APERO_BIN_PATH/$INSTALL_SCRIPT_REL"

# -----------------------------------------------------------------------------
# Check if the setup script file actually exists
# -----------------------------------------------------------------------------
if [ ! -f "$SETUP_SCRIPT" ]; then
  echo ""
  echo "Error: Setup script for instrument '$INSTRUMENT_NAME' not found at:"
  echo "  $SETUP_SCRIPT"
  echo ""
  echo "The instrument is defined in $INSTRUMENT_FILE but the script file is missing."
  echo ""
  echo "Available instruments in $INSTRUMENT_FILE:"
  grep -E '^\[[^]]+\]' "$INSTRUMENT_FILE" | sed 's/^\[//; s/\]$//'
  echo ""
  exit 1
fi

read -r -d '' SNIPPET <<EOF
# Source $INSTRUMENT_NAME profile if present
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
    echo "Installing instrument profile for '$INSTRUMENT_NAME' into $TARGET"
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

