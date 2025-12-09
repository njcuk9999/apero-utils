#!/bin/bash

# Set the project ID (it may change in future)
APERO_PROJECT_ID="6102120"

# Set the apero bin path
APERO_BIN_PATH="/project/$APERO_PROJECT_ID/apero/apero_bin"

# Set the path to instruments ini
INSTRUMENT_FILE="$APERO_BIN_PATH/apero_instruments.ini"

# Track the last valid git path
LAST_VALID_GIT_PATH=""

echo "Reading instruments from: $INSTRUMENT_FILE"
echo

# Loop through instruments in key=value format
while IFS='=' read -r INSTRUMENT SCRIPT; do
    # Skip empty or comment lines
    [[ -z "$INSTRUMENT" || "$INSTRUMENT" =~ ^# ]] && continue

    INSTRUMENT_BIN="${INSTRUMENT}_bin"

    GIT_PATH="/project/$APERO_PROJECT_ID/apero/${INSTRUMENT_BIN}/scripts/apero-utils/alliance/apero_bin"

    echo "Checking: $GIT_PATH"

    if [[ -d "$GIT_PATH" ]]; then
        echo "→ Found. Entering directory..."
        cd "$GIT_PATH" || continue
        git pull
        LAST_VALID_GIT_PATH="$GIT_PATH"
        echo "✓ Updated $INSTRUMENT"
    else
        echo "✗ Path does not exist. Skipping $INSTRUMENT"
    fi

    echo
done < "$INSTRUMENT_FILE"

echo "Done looping over instruments."
echo

# After loop: ensure we have a valid git path
if [[ -z "$LAST_VALID_GIT_PATH" ]]; then
    echo "ERROR: No valid GIT_PATH directories found. Nothing to copy."
    exit 1
fi

echo "Last valid GIT_PATH: $LAST_VALID_GIT_PATH"
echo "Copying updated files back to APERO_BIN_PATH..."

cd "$APERO_BIN_PATH" || exit 1

cp -v "$LAST_VALID_GIT_PATH"/* "$APERO_BIN_PATH"

echo "Sync complete."

