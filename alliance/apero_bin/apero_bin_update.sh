#!/bin/bash

# Get core source script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $SCRIPT_DIR/apero_core.sh
# Set the path to instruments conf
INSTRUMENT_FILE="$APERO_BIN_PATH/apero_instruments.conf"

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
    echo "  ./apero_bin_update.sh"
    echo ""
    echo "Notes:"
    echo "  - Ensure APERO_PROJECT_ID is set correctly."
    echo "  - The script expects $INSTRUMENT_FILE sections [instrument] and an install_script line."
    echo "  - Empty lines and lines starting with # are ignored."
    echo "  - Requires git to be installed and accessible in PATH."
    echo ""
    echo "Example:"
    echo "  ./apero_bin_update.sh"
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

if [[ ! -f "$INSTRUMENT_FILE" ]]; then
    echo "ERROR: Instruments config not found: $INSTRUMENT_FILE"
    exit 1
fi

echo "Reading instruments from: $INSTRUMENT_FILE"
echo

CURRENT_INSTRUMENT=""
INSTALL_SCRIPT=""

process_instrument() {
    local instrument="$1"
    local install_script="$2"
    if [[ -z "$instrument" ]]; then
        return
    fi

    local INSTRUMENT_BIN="${instrument}_bin"
    local GIT_PATH="/project/$APERO_PROJECT_ID/apero/${INSTRUMENT_BIN}/scripts/apero-utils/alliance/apero_bin"

    echo "Instrument: $instrument"
    echo "Install script: $install_script"
    echo "Checking: $GIT_PATH"

    if [[ -d "$GIT_PATH" ]]; then
        echo "→ Found. Entering directory..."
        cd "$GIT_PATH" || return
        git pull || { echo "✗ git pull failed"; echo; return; }
        LAST_VALID_GIT_PATH="$GIT_PATH"
        echo "✓ Updated $instrument"
    else
        echo "✗ Path does not exist. Skipping $instrument"
    fi

    echo
}

# Parse the INSTRUMENT_FILE
while IFS= read -r line || [ -n "$line" ]; do
    line="$(echo "$line" | sed 's/^\s\+//;s/\s\+$//')"
    [[ -z "$line" || "$line" =~ ^# ]] && continue

    if [[ "$line" =~ ^\[(.+)\]$ ]]; then
        # new instrument section
        if [[ -n "$CURRENT_INSTRUMENT" ]]; then
            process_instrument "$CURRENT_INSTRUMENT" "$INSTALL_SCRIPT"
        fi
        CURRENT_INSTRUMENT="${BASH_REMATCH[1]}"
        INSTALL_SCRIPT=""
        echo "----------------------------"
        echo "Section: $CURRENT_INSTRUMENT"
        continue
    fi

    if [[ "$line" =~ ^install_script=(.+)$ ]]; then
        INSTALL_SCRIPT="${BASH_REMATCH[1]}"
        continue
    fi

done < "$INSTRUMENT_FILE"

# process last instrument
if [[ -n "$CURRENT_INSTRUMENT" ]]; then
    process_instrument "$CURRENT_INSTRUMENT" "$INSTALL_SCRIPT"
fi

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

# -----------------------------------------------------------------------------
# Copy docs from apero-utils/alliance/docs to $APERO_PATH/docs/
# -----------------------------------------------------------------------------
# Source docs path (assume apero-utils is under $APERO_PATH)
DOCS_SRC="$APERO_PATH/apero-utils/alliance/docs"
DOCS_DEST_DIR="$APERO_PATH/docs"

if [[ -d "$DOCS_SRC" ]]; then
    echo "Copying docs from $DOCS_SRC to $DOCS_DEST_DIR"
    mkdir -p "$DOCS_DEST_DIR"
    # Copy contents (preserve attributes), avoid nesting by copying contents
    cp -a "$DOCS_SRC/." "$DOCS_DEST_DIR/" || {
        echo "Warning: failed to copy docs from $DOCS_SRC to $DOCS_DEST_DIR"
    }
    echo "Docs copied."
else
    echo "Note: docs not found at $DOCS_SRC - skipping docs copy."
fi

echo "Sync complete."

