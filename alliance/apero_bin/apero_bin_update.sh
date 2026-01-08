#!/bin/bash

# Get core source script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $SCRIPT_DIR/apero_core.sh
# Set the path to instruments ini
INSTRUMENT_FILE="$APERO_BIN_PATH/apero_instruments.ini"

# -----------------------------------------------------------------------------
# Function: Show help message
# -----------------------------------------------------------------------------
show_help() {
    echo ""
    echo "APERO Instrument Git Updater"
    echo "----------------------------"
    echo "This script loops through instruments defined in:"
    echo "  $INSTRUMENT_FILE"
    echo "and attempts to update the corresponding git directories."
    echo ""
    echo "It then copies the last successfully updated instrument files"
    echo "back to the APERO bin directory:"
    echo "  $APERO_BIN_PATH"
    echo ""
    echo "Usage:"
    echo "  ./this_script.sh"
    echo ""
    echo "Notes:"
    echo "  - Ensure APERO_PROJECT_ID is set correctly."
    echo "  - The script expects each line of $INSTRUMENT_FILE to be in key=value format:"
    echo "      instrument_name=script_name"
    echo "  - Empty lines and lines starting with # are ignored."
    echo "  - Requires git to be installed and accessible in PATH."
    echo ""
    echo "Example:"
    echo "  ./apero_update_instruments.sh"
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

cp "$LAST_VALID_GIT_PATH"/* "$APERO_BIN_PATH"

# Set permissions: 750 for .sh files, 640 for config files
chmod 750 "$APERO_BIN_PATH"/*.sh
chmod 640 "$APERO_BIN_PATH"/*.conf "$APERO_BIN_PATH"/*.ini 2>/dev/null || true

echo "Sync complete."

