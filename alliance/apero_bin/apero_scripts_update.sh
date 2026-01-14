#!/bin/bash

# -----------------------------------------------------------------------------
# APERO Scripts Git Updater
# -----------------------------------------------------------------------------

# Get core source script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/apero_core.sh"

# Instruments configuration (bash-friendly .conf)
INSTRUMENT_FILE="$APERO_BIN_PATH/apero_instruments.conf"

# DEBUG: set DEBUG=1 in environment to see parsing debug

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
    echo "  DEBUG=1 ./apero_scripts_update.sh    # show parsing debug"
    echo ""
    echo "Notes:"
    echo "  - Ensure APERO_PROJECT_ID is set correctly (apero_core.sh)."
    echo "  - The config format is bash-friendly, sections [instrument],"
    echo "    with install_script and multiple repo lines:"
    echo "      repo=dir_name branch=branch_name"
    echo "  - Requires git to be installed and accessible in PATH."
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

# Helper: debug print if DEBUG=1
dbg() {
    if [[ "${DEBUG:-0}" -eq 1 ]]; then
        echo "DEBUG: $*"
    fi
}

process_instrument() {
    local instrument="$1"
    local install_script="$2"
    shift 2
    local repos=("$@")

    if [[ -z "$instrument" ]]; then
        return
    fi

    echo "=================================================="
    echo "Instrument: $instrument"
    echo "=================================================="

    local INSTRUMENT_BIN="${instrument}_bin"
    # Prefer APERO_PATH (set by apero_core.sh). Fallback to constructed /project path
    local BASE_PATH="${APERO_PATH:-/project/$APERO_PROJECT_ID/apero}"
    local SCRIPT_PATH="${BASE_PATH}/${INSTRUMENT_BIN}/scripts"

    # If scripts path doesn't exist, try the bin root as fallback
    if [[ ! -d "$SCRIPT_PATH" && -d "${BASE_PATH}/${INSTRUMENT_BIN}" ]]; then
        dbg "note: scripts subdir not found; using ${BASE_PATH}/${INSTRUMENT_BIN} as script path"
        SCRIPT_PATH="${BASE_PATH}/${INSTRUMENT_BIN}"
    fi

    # If still not found, try to discover candidate directories under BASE_PATH
    if [[ ! -d "$SCRIPT_PATH" ]]; then
        dbg "attempting to discover instrument directory under $BASE_PATH"
        # look for directories matching instrument name, prefer ones ending with _bin
        candidate=$(find "$BASE_PATH" -maxdepth 3 -type d -iname "*${instrument}*bin" -print -quit 2>/dev/null || true)
        if [[ -z "$candidate" ]]; then
            candidate=$(find "$BASE_PATH" -maxdepth 4 -type d -iname "*${instrument}*" -print -quit 2>/dev/null || true)
        fi
        if [[ -n "$candidate" ]]; then
            dbg "Found candidate instrument dir: $candidate"
            if [[ -d "$candidate/scripts" ]]; then
                SCRIPT_PATH="$candidate/scripts"
            else
                SCRIPT_PATH="$candidate"
            fi
        fi
    fi

    if [[ ! -d "$SCRIPT_PATH" ]]; then
        echo "✗ Scripts path does not exist: $SCRIPT_PATH"
        echo
        return
    fi

    echo "Scripts path: $SCRIPT_PATH"
    echo

    cd "$SCRIPT_PATH" || {
        echo "ERROR: Cannot cd to $SCRIPT_PATH"
        echo
        return
    }

    if [[ ${#repos[@]} -eq 0 ]]; then
        echo "  (No repos configured for instrument: $instrument)"
        echo
        return
    fi

    for entry in "${repos[@]}"; do
        # entry format: dir_name|branch_name
        local repo_dir="${entry%%|*}"
        local branch_name="${entry##*|}"
        local REPO_PATH="$SCRIPT_PATH/$repo_dir"

        echo "--------------------------------------------------"
        echo "Repo: $repo_dir (branch: $branch_name)"
        echo "--------------------------------------------------"

        if [[ ! -d "$REPO_PATH/.git" ]]; then
            echo "✗ Not a git repo: $REPO_PATH"
            echo
            continue
        fi

        echo "Fetching..."
        (cd "$REPO_PATH" && git fetch --all --prune) || { echo "✗ git fetch failed"; echo; continue; }

        echo "Checking out branch: $branch_name"
        (cd "$REPO_PATH" && git checkout "$branch_name") || { echo "✗ git checkout failed"; echo; continue; }

        echo "Pulling..."
        (cd "$REPO_PATH" && git pull --ff-only) || { echo "✗ git pull failed"; echo; continue; }

        echo "✓ Updated $repo_dir"
        echo
    done

    echo
}

# Parse the INSTRUMENT_FILE
while IFS= read -r line || [ -n "$line" ]; do
    # trim whitespace (POSIX safe)
    line="$(printf '%s' "$line" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')"
    dbg "parsed line='${line}'"
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
        continue
    fi

    if [[ "$line" =~ ^install_script=(.+)$ ]]; then
        INSTALL_SCRIPT="${BASH_REMATCH[1]}"
        continue
    fi

    # robust parsing for repo lines: repo=<name> [other tokens] branch=<branch>
    if [[ "$line" == repo=* ]]; then
        rest="${line#repo=}"
        # first token is repo dir
        repo_dir="$(printf '%s' "$rest" | awk '{print $1}')"
        # extract branch=NAME if present
        branch_name="$(printf '%s' "$rest" | sed -n 's/.*branch=\([^[:space:]]*\).*/\1/p')"
        # default branch if none provided
        if [[ -z "$branch_name" ]]; then
            branch_name="master"
        fi
        dbg "Parsed repo: $repo_dir  branch: $branch_name"
        REPOS+=("${repo_dir}|${branch_name}")
        continue
    fi

done < "$INSTRUMENT_FILE"

# process last instrument
if [[ -n "$CURRENT_INSTRUMENT" ]]; then
    process_instrument "$CURRENT_INSTRUMENT" "$INSTALL_SCRIPT" "${REPOS[@]}"
fi

echo "Sync complete."

