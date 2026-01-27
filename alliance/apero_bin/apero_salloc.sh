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
  -p, --prompt, --run   Prompt the user for confirmation before running salloc
  -h, --help            Show this help message and exit

If an option is not provided on the command line the script will prompt for it
interactively (except the final confirmation which is only asked when
--prompt is provided). Without --prompt the script will not run salloc and
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
        -p|--prompt|--run)
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
# Default values (used only if user accepts empty input at prompts)
# Do NOT assign these defaults to the variables here — only prompt if not
# supplied on the command line (unless --prompt was provided, see below).
# -----------------------------------------------------------------------------
DEFAULT_TIME="4:00:00"
DEFAULT_NODES="1"
DEFAULT_CPUS="1"
DEFAULT_MEM="4096M"

# -----------------------------------------------------------------------------
# Build account list from conf file (we need it both for interactive and
# non-interactive/--prompt modes)
# -----------------------------------------------------------------------------

echo
echo "Available $APERO_SERVER accounts:"
echo

i=1
declare -a ACCOUNT_LIST
current_group=""
while IFS= read -r line; do
    if [[ "$line" =~ ^\[(.+)\]$ ]]; then
        current_group="${BASH_REMATCH[1]}"
        continue
    fi
    [[ -z "$line" ]] && continue
    ACCOUNT_LIST[$i]="$line"
    echo "  $i) $line   (group: $current_group)"
    ((i++))

done < "$CONF_FILE"

if (( i == 1 )); then
    echo "No accounts found!"
    exit 1
fi

# -----------------------------------------------------------------------------
# Two modes:
#  - PROMPT=1: user requested confirmation of running the salloc command.
#              In this mode we DO NOT interactively ask for missing options;
#              instead we auto-fill them with defaults (and pick first
#              account if none provided), then ask the single confirmation
#              "Run salloc [Y]es or [N]o".
#  - PROMPT=0: interactive mode for building the command: ask for any
#              missing options (unless provided on CLI). After collecting
#              values we display the command and exit (do not run it).
# -----------------------------------------------------------------------------

if [[ $PROMPT -eq 1 ]]; then
    # Non-interactive option collection: fill missing values with defaults
    TIME=${TIME:-$DEFAULT_TIME}
    NODES=${NODES:-$DEFAULT_NODES}
    CPUS=${CPUS:-$DEFAULT_CPUS}
    MEM=${MEM:-$DEFAULT_MEM}

    # If account not provided, pick the first account and inform the user
    if [[ -z "$ACCOUNT" ]]; then
        ACCOUNT="${ACCOUNT_LIST[1]}"
        echo "No account provided; defaulting to: $ACCOUNT"
    else
        echo "Using account: $ACCOUNT"
    fi

    # X11 flag handling when non-interactive: respect WANT_X11 if set, else none
    if [[ "$WANT_X11" == "Y" ]]; then
        X11_FLAG="--x11"
    else
        X11_FLAG=""
    fi

else
    # Interactive mode: prompt for any missing options
    if [[ -z "${TIME+x}" || -z "$TIME" ]]; then
        read -p "Enter time (HH:MM:SS) [${DEFAULT_TIME}]: " TIME
        TIME=${TIME:-$DEFAULT_TIME}
    else
        echo "Using time: $TIME"
    fi

    if [[ -z "${NODES+x}" || -z "$NODES" ]]; then
        read -p "Enter number of nodes [${DEFAULT_NODES}]: " NODES
        NODES=${NODES:-$DEFAULT_NODES}
    else
        echo "Using nodes: $NODES"
    fi

    if [[ -z "${CPUS+x}" || -z "$CPUS" ]]; then
        read -p "Enter number of CPUs per task [${DEFAULT_CPUS}]: " CPUS
        CPUS=${CPUS:-$DEFAULT_CPUS}
    else
        echo "Using cpus: $CPUS"
    fi

    if [[ -z "${MEM+x}" || -z "$MEM" ]]; then
        read -p "Enter mem per CPU (e.g., 4096M) [${DEFAULT_MEM}]: " MEM
        MEM=${MEM:-$DEFAULT_MEM}
    else
        echo "Using mem: $MEM"
    fi

    # Interactive account selection if not provided
    if [[ -n "$ACCOUNT" ]]; then
        echo "Using account: $ACCOUNT"
    else
        echo
        read -p "Select account by number: " CHOICE
        ACCOUNT="${ACCOUNT_LIST[$CHOICE]}"
        if [[ -z "$ACCOUNT" ]]; then
            echo "Invalid selection."
            exit 1
        fi
    fi

    # X11 interactive prompt
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

# Final action: if PROMPT=1 ask to run; if PROMPT=0 do not run (interactive
# mode already collected options and displayed the command).
if [[ $PROMPT -eq 1 ]]; then
    read -p "Run salloc [Y]es or [N]o: " CONFIRM
    CONFIRM=${CONFIRM:-N}
    if [[ ! "$CONFIRM" =~ ^[Yy] ]]; then
        echo "Aborted by user."
        exit 1
    fi
    eval "$COMMAND"
else
    echo "No --prompt flag provided; command shown but not executed."
    echo "If you want to confirm and run the command add --prompt to the command line."
    exit 0
fi
