#!/bin/bash

# -----------------------------------------------------------------------------
# APERO Scripts Git Updater
# -----------------------------------------------------------------------------

# Get core source script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/apero_core.sh"

# Instruments configuration (bash-friendly .conf)
INSTRUMENT_FILE="$APERO_BIN_PATH/apero_instruments.conf"

# -----------------------------------------------------------------------------
# Function: Show help message
# -----------------------------------------------------------------------------
show_help() {
    echo ""
    echo "APERO Scripts Git Updater"
    echo "-------------------------"
    echo "Reads instruments and repositories from:"
    echo "  $INSTRUMENT_FILE"
    echo "Fetches, checks out branches, and pulls git repositories for each instrument."
    echo ""
    echo "Usage:"
    echo "  ./apero_scripts_update.sh"
    echo ""
    echo "Notes:"
    echo "  - Ensure APERO_PROJECT_ID is set correctly."
    echo "  - The config format is bash-friendly, sections [instrument],"
    echo "    with install_script and multiple repo lines:"
    echo "      repo=dir_name branch=branch_name"
    echo "  - Requires git to be installed and accessible in PATH."
    echo ""
    echo "Example:"
    echo "  ./apero_scripts_update.sh"
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

if [[ ! -f "$INSTRUMENT_FILE" ]]; then
    echo "ERROR: Instruments config not found: $INSTRUMENT_FILE"
    exit 1
fi

echo "Reading instruments from: $INSTRUMENT_FILE"
echo

CURRENT_INSTRUMENT=""
INSTALL_SCRIPT=""
REPOS=()

process_instrument() {
    local instrument="$1"
    local install_script="$2"
    shift 2
    local repos=("$@")

    if [[ -z "$instrument" ]]; then
        return
    fi

    echo "Instrument: $instrument"
    echo "Install script: $install_script"

    local INSTRUMENT_BIN="${instrument}_bin"
    local SCRIPT_PATH="/project/$APERO_PROJECT_ID/apero/${INSTRUMENT_BIN}/scripts"

    if [[ ! -d "$SCRIPT_PATH" ]]; then
        echo "✗ Scripts path does not exist: $SCRIPT_PATH"
        echo
        return
    fi

    echo "→ Entering scripts path: $SCRIPT_PATH"
    cd "$SCRIPT_PATH" || {
        echo "ERROR: Cannot cd to $SCRIPT_PATH"
        echo
        return
    }

    for entry in "${repos[@]}"; do
        # entry format: dir_name|branch_name
        local repo_dir="${entry%%|*}"
        local branch_name="${entry##*|}"
        local REPO_PATH="$SCRIPT_PATH/$repo_dir"

        echo "  Repo: $repo_dir (branch: $branch_name)"
        if [[ ! -d "$REPO_PATH/.git" ]]; then
            echo "    ✗ Not a git repo: $REPO_PATH"
            continue
        fi

        echo "    → Fetching..."
        (cd "$REPO_PATH" && git fetch) || { echo "    ✗ git fetch failed"; continue; }
        echo "    → Checking out branch: $branch_name"
        (cd "$REPO_PATH" && git checkout "$branch_name") || { echo "    ✗ git checkout failed"; continue; }
        echo "    → Pulling..."
        (cd "$REPO_PATH" && git pull) || { echo "    ✗ git pull failed"; continue; }
        echo "    ✓ Updated $repo_dir"
    done

    echo
}

# Parse the INSTRUMENT_FILE
while IFS= read -r line || [ -n "$line" ]; do
    # trim whitespace
    line="$(echo "$line" | sed 's/^\s\+//;s/\s\+$//')"
    # skip comments and empty lines
    [[ -z "$line" || "$line" =~ ^# ]] && continue

    if [[ "$line" =~ ^\[(.+)\]$ ]]; then
        # new instrument section
        # process previous instrument
        if [[ -n "$CURRENT_INSTRUMENT" ]]; then
            process_instrument "$CURRENT_INSTRUMENT" "$INSTALL_SCRIPT" "${REPOS[@]}"
        fi
        CURRENT_INSTRUMENT="${BASH_REMATCH[1]}"
        INSTALL_SCRIPT=""
        REPOS=()
        echo "----------------------------"
        echo "Section: $CURRENT_INSTRUMENT"
        continue
    fi

    if [[ "$line" =~ ^install_script=(.+)$ ]]; then
        INSTALL_SCRIPT="${BASH_REMATCH[1]}"
        continue
    fi

    if [[ "$line" =~ ^repo=([^[:space:]]+)\s+branch=([^[:space:]]+)$ ]]; then
        repo_dir="${BASH_REMATCH[1]}"
        branch_name="${BASH_REMATCH[2]}"
        REPOS+=("${repo_dir}|${branch_name}")
        continue
    fi

done < "$INSTRUMENT_FILE"

# process last instrument
if [[ -n "$CURRENT_INSTRUMENT" ]]; then
    process_instrument "$CURRENT_INSTRUMENT" "$INSTALL_SCRIPT" "${REPOS[@]}"
fi

echo "Sync complete."

