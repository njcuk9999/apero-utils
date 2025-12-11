#!/bin/bash

# Get core source script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $SCRIPT_DIR/apero_core.sh

# define conf file
CONF_FILE="$APERO_BIN_PATH/apero_groups.conf"

show_help() {
    cat << EOF

Usage: apero_salloc.sh

This script helps you launch a SLURM salloc session interactively.

It will ask for the following options:

  - Time          : Job time in HH:MM:SS format (default: 4:00:00)
  - CPUs          : Number of CPUs per task (default: 20)
  - Nodes         : Number of nodes (default: 1)
  - Memory        : Memory per CPU, e.g., 4096M (default: 4096M)
  - Account       : Select from accounts defined in your apero_users.conf
  - X11           : Optional interactive X11 session
  - Confirmation  : Asks before running the salloc command

Your accounts are read from the conf file:
    $CONF_FILE

Environment variables used:
  - APERO_SERVER     : default server/account (if used)

Example usage:

    apero_salloc.sh

Pressing Ctrl+C at any prompt will abort the script.

EOF
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
# Test if conf file exists
# -----------------------------------------------------------------------------
if [[ ! -f "$CONF_FILE" ]]; then
    echo "Error: $CONF_FILE not found!"
    exit 1
fi

echo "=== SLURM salloc launcher ==="

# -----------------------------------------------------------------------------
# Ask about time allocation
# -----------------------------------------------------------------------------
read -p "Enter time (HH:MM:SS) [default 4:00:00]: " TIME
TIME=${TIME:-4:00:00}

# -----------------------------------------------------------------------------
# Ask about Nodes
# -----------------------------------------------------------------------------
read -p "Enter number of nodes [default 1]: " NODES
NODES=${NODES:-1}

# -----------------------------------------------------------------------------
# Ask about CPUS
# -----------------------------------------------------------------------------
read -p "Enter number of CPUs per task [default 1]: " CPUS
CPUS=${CPUS:-1}

# -----------------------------------------------------------------------------
# Ask about Memory
# -----------------------------------------------------------------------------
read -p "Enter mem per CPU (e.g., 4096M) [default 4096M]: " MEM
MEM=${MEM:-4096M}

# -----------------------------------------------------------------------------
# Ask about user account
# -----------------------------------------------------------------------------
echo
echo "Available $APERO_SERVER accounts:"
echo

# Parse groups and map them to env variables
i=1
declare -a ACCOUNT_LIST

current_group=""

while IFS= read -r line; do
    # Detect group header
    if [[ "$line" =~ ^\[(.+)\]$ ]]; then
        current_group="${BASH_REMATCH[1]}"
        continue
    fi

    # Skip blank lines
    [[ -z "$line" ]] && continue

    # Treat every non-header line as an account entry
    ACCOUNT_LIST[$i]="$line"
    echo "  $i) $line   (group: $current_group)"
    ((i++))

done < "$CONF_FILE"

if (( i == 1 )); then
    echo "No accounts found!"
    exit 1
fi

echo
read -p "Select account by number: " CHOICE

ACCOUNT="${ACCOUNT_LIST[$CHOICE]}"

if [[ -z "$ACCOUNT" ]]; then
    echo "Invalid selection."
    exit 1
fi

# -----------------------------------------------------------------------------
# Ask about interactive session
# -----------------------------------------------------------------------------
# Ask whether user wants an interactive X11 session
X11_FLAG=""

echo
read -p "Request an interactive X11 session? [Y/N]: " WANT_X11
WANT_X11=${WANT_X11:-N}

if [[ "$WANT_X11" =~ ^[Yy]$ ]]; then
    X11_FLAG="--x11"
fi

# -----------------------------------------------------------------------------
# Make salloc command
# -----------------------------------------------------------------------------

COMMAND="salloc --time=$TIME --cpus-per-task=$CPUS --nodes=$NODES --mem-per-cpu=$MEM --account=$ACCOUNT $X11_FLAG"

echo
echo "The following salloc command will be run:"
echo
echo ">> $COMMAND"
echo

read -p "Run this command? [Y/n]: " CONFIRM
CONFIRM=${CONFIRM:-Y}

if [[ ! "$CONFIRM" =~ ^[Yy]$ ]]; then
    echo "Aborted by user."
    exit 1
fi

eval "$COMMAND"