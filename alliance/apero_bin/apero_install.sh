#!/bin/bash


# Note this code is only to be run once per user
#    Do not add anything here to change setups
#    This just puts stuff in the ~/.bashrc
#    Please look at apero_instruement.ini to see list of setup scripts

# Set the project ID (it may change in future)
APERO_PROJECT_ID="6102120"

# Find the bash file to push this into
TARGET="$HOME/.bashrc"

# Set the apero bin path
APERO_BIN_PATH="/project/$APERO_PROJECT_ID/apero/apero_bin"

# Set the path to profiles.ini
INSTRUMENT_FILE="$APERO_BIN_PATH/apero_instruments.ini"


# Check if apero_instruments.ini exists
if [ ! -f "$INSTRUMENT_FILE" ]; then
  echo "Error: $INSTRUMENT_FILE file not found."
  return 1
fi

# Check if instrument name is provided as an argument
if [ -z "$1" ]; then
  echo ""
  echo "Usage: apero_install.sh <instrument_name>"
  echo ""
  echo "Available instruments are:"
  grep -o '^[^=]*' $INSTRUMENT_FILE
  echo ""
  return 1
fi

# print that we found it (may remove later)
if [ "$2" = "--debug" ]; then
  echo "Found: INSTRUMENT_FILE=$INSTRUMENT_FILE"
fi

# Read profiles.ini and find the path for the provided profile
INSTRUMENT_PATH=$(grep "^$1=" $INSTRUMENT_FILE | cut -d'=' -f2)

# Check if profile path exists
if [ -z "$INSTRUMENT_FILE" ]; then
  echo ""
  echo "Profile $1 not found in $INSTRUMENT_FILE."
  echo ""
  echo "Available profiles are:"
  grep -o '^[^=]*' $INSTRUMENT_FILE
  echo ""
  echo "Or run apero_setup.py to create a new profile"
  echo ""
  return 1
fi

# print that we found it (may remove later)
if [ "$2" = "--debug" ]; then
  echo "Found: INSTRUMENT_PATH=$INSTRUMENT_PATH"
fi

# Determine the setup script to run based on OS
SETUP_SCRIPT="$APERO_BIN_PATH/$INSTRUMENT_PATH"

read -r -d '' SNIPPET <<EOF
# Source instrument profile if present
if [ -f "$SETUP_SCRIPT" ]; then
    source "$SETUP_SCRIPT"
fi
EOF

# print that we found it (may remove later)
if [ "$2" = "--debug" ]; then
  echo "SNIPPET = $SNIPPET"
fi

# Append only if not present
if ! grep -Fq "$SETUP_SCRIPT" "$TARGET"; then
    printf "\n%s\n" "$SNIPPET" >> "$TARGET"
    source "$SETUP_SCRIPT"
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




