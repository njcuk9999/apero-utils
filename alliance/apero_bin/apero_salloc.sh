#!/bin/bash

# Get core source script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $SCRIPT_DIR/apero_core.sh

# define conf file
CONF_FILE="$APERO_BIN_PATH/apero_groups.conf"

show_help() {
    cat << EOF

Usage: apero_salloc.sh [options]

This script helps you launch a SLURM salloc session interactively or via
command-line options.

Options:
  -t, --time TIME       Job time in HH:MM:SS format (default: 4:00:00)
  -c, --cpus CPUS       Number of CPUs per task (default: 1)
  -n, --nodes NODES     Number of nodes (default: 1)
  -m, --mem MEM         Memory per CPU, e.g., 4096M (default: 4096M)
  -a, --account ACCOUNT Select account by name from the conf file
      --x11            Request an interactive X11 session
      --prompt         Prompt the user for confirmation before running salloc
  -h, --help            Show this help message and exit

If an option is not provided on the command line the script will prompt for it
interactively (except the final confirmation which is only asked when
`--prompt` is provided). Without `--prompt` the script will not run salloc and
will exit after showing the command to be run.

Your accounts are read from the conf file:
    $CONF_FILE

Environment variables used:
  - APERO_SERVER     : default server/account (if used)

Example usage:

    apero_salloc.sh --time 04:00:00 --cpus 20 --nodes 1 --mem 4096M --account myacct --prompt

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
# Parse command line arguments (short and long)
# -----------------------------------------------------------------------------
PROMPT=0
WANT_X11="N"
# Variables left empty will be prompted for later
# TIME, NODES, CPUS, MEM, ACCOUNT may be set here
while [[ $# -gt 0 ]]; do
    case "$1" in
        -t|--time)
            TIME="$2"; shift 2;;
        -c|--cpus)
            CPUS="$2"; shift 2;;
        -n|--nodes)
            NODES="$2"; shift 2;;
        -m|--mem)
            MEM="$2"; shift 2;;
        -a|--account)
            ACCOUNT="$2"; shift 2;;
        --x11)
            WANT_X11="Y"; shift;;
        --prompt)
            PROMPT=1; shift;;
        -h|--help)
            show_help; exit 0;;
        *)
            echo "*************************"
            echo "ERROR: Unknown option: $1"
            echo "*************************"
            show_help
            exit 1;;
    esac
done

# -----------------------------------------------------------------------------
# Test if conf file exists
# -----------------------------------------------------------------------------
if [[ ! -f "$CONF_FILE" ]]; then
    echo "*************************"
    echo "ERROR: $CONF_FILE not found!"
    echo "*************************"
    exit 1
fi

echo "=== SLURM salloc launcher ==="

# -----------------------------------------------------------------------------
# Set defaults for values not provided on the command line
# -----------------------------------------------------------------------------
TIME=${TIME:-4:00:00}
NODES=${NODES:-1}
CPUS=${CPUS:-1}
MEM=${MEM:-4096M}

# -----------------------------------------------------------------------------
# Ask about time allocation if not provided
# -----------------------------------------------------------------------------
if [[ -z "$TIME" ]]; then
    read -p "Enter time (HH:MM:SS) [default 4:00:00]: " TIME
    TIME=${TIME:-4:00:00}
else
    echo "Using time: $TIME"
fi

# -----------------------------------------------------------------------------
# Ask about Nodes if not provided
# -----------------------------------------------------------------------------
if [[ -z "$NODES" ]]; then
    read -p "Enter number of nodes [default 1]: " NODES
    NODES=${NODES:-1}
else
    echo "Using nodes: $NODES"
fi

# -----------------------------------------------------------------------------
# Ask about CPUS if not provided
# -----------------------------------------------------------------------------
if [[ -z "$CPUS" ]]; then
    read -p "Enter number of CPUs per task [default 1]: " CPUS
    CPUS=${CPUS:-1}
else
    echo "Using cpus: $CPUS"
fi

# -----------------------------------------------------------------------------
# Ask about Memory if not provided
# -----------------------------------------------------------------------------
if [[ -z "$MEM" ]]; then
    read -p "Enter mem per CPU (e.g., 4096M) [default 4096M]: " MEM
    MEM=${MEM:-4096M}
else
    echo "Using mem: $MEM"
fi

# -----------------------------------------------------------------------------
# Ask about user account - build ACCOUNT_LIST from conf file
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

# If account provided via CLI, validate it exists in the ACCOUNT_LIST values
if [[ -n "$ACCOUNT" ]]; then
    found=0
    for ((j=1;j<i;j++)); do
        if [[ "${ACCOUNT_LIST[$j]}" == "$ACCOUNT" ]]; then
            found=1
            break
        fi
    done
    if [[ $found -ne 1 ]]; then
        echo "*************************"
        echo "ERROR: Account '$ACCOUNT' not found in $CONF_FILE"
        echo "*************************"
        exit 1
    fi
else
    # interactive selection
    echo
    read -p "Select account by number: " CHOICE
    ACCOUNT="${ACCOUNT_LIST[$CHOICE]}"

    if [[ -z "$ACCOUNT" ]]; then
        echo "Invalid selection."
        exit 1
    fi
fi

# -----------------------------------------------------------------------------
# Ask about interactive session
# -----------------------------------------------------------------------------

if [[ "$WANT_X11" == "Y" ]]; then
    X11_FLAG="--x11"
else
    echo
    read -p "Request an interactive X11 session? [Y/N]: " WANT_X11
    WANT_X11=${WANT_X11:-N}

    if [[ "$WANT_X11" =~ ^[Yy]$ ]]; then
        X11_FLAG="--x11"
    else
        X11_FLAG=""
    fi
fi

# -----------------------------------------------------------------------------
# Make salloc command
# -----------------------------------------------------------------------------

COMMAND="salloc --time=$TIME --cpus-per-task=$CPUS --nodes=$NODES --mem-per-cpu=$MEM --account=$ACCOUNT $X11_FLAG"

echo
echo "The following salloc command will be run:"
echo
echo "=================================================="
echo ">> $COMMAND"
echo "=================================================="
echo

# If --prompt was provided ask the user, otherwise do not run and inform the user
if [[ $PROMPT -eq 1 ]]; then
    read -p "Run this command? [Y/n]: " CONFIRM
    CONFIRM=${CONFIRM:-Y}

    if [[ ! "$CONFIRM" =~ ^[Yy]$ ]]; then
        echo "Aborted by user."
        exit 1
    fi

    eval "$COMMAND"
else
    echo "No --prompt flag provided; not running salloc."
    echo "If you want to be prompted to run the command add --prompt to the command line."
    exit 0
fi

