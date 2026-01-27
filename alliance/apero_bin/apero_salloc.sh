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

# If --prompt was provided ask the user immediately whether to run; this
# appears before any other interactive prompts per user request.
RUN_CONFIRMED=0
if [[ $PROMPT -eq 1 ]]; then
    if [[ -c /dev/tty ]]; then
        read -p "Run salloc [Y]es or [N]o: " FIRST_CONFIRM </dev/tty
    else
        read -p "Run salloc [Y]es or [N]o: " FIRST_CONFIRM
    fi
    FIRST_CONFIRM=${FIRST_CONFIRM:-N}
    if [[ ! "$FIRST_CONFIRM" =~ ^[Yy] ]]; then
        echo "Aborted by user."
        exit 1
    fi
    RUN_CONFIRMED=1
fi

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

i=1
declare -a ACCOUNT_LIST
declare -a ACCOUNT_GROUP
current_group=""
while IFS= read -r line; do
    if [[ "$line" =~ ^\[(.+)\]$ ]]; then
        current_group="${BASH_REMATCH[1]}"
        continue
    fi
    [[ -z "$line" ]] && continue
    ACCOUNT_LIST[$i]="$line"
    ACCOUNT_GROUP[$i]="$current_group"
    ((i++))
done < "$CONF_FILE"

if (( i == 1 )); then
    echo "No accounts found!"
    exit 1
fi

# -----------------------------------------------------------------------------
# helper to prompt only when interactive (stdin or stdout is a tty or /dev/tty exists)
is_interactive=0
if [ -t 0 ] || [ -t 1 ] || [ -t 2 ] || [ -c /dev/tty ]; then
    is_interactive=1
fi

prompt_default() {
    # args: varname prompt_text default
    local __varname="$1"; shift
    local __prompt="$1"; shift
    local __default="$1"; shift

    if [[ $is_interactive -eq 1 ]]; then
        # Use /dev/tty when available so prompts show when stdin is redirected
        if [[ -c /dev/tty ]]; then
            read -p "${__prompt} [${__default}]: " __input </dev/tty
        else
            read -p "${__prompt} [${__default}]: " __input
        fi
        if [[ -z "$__input" ]]; then
            eval "${__varname}=\"${__default}\""
        else
            # assign the input
            eval "${__varname}=\"$__input\""
        fi
    else
        # Non-interactive: auto-fill default
        eval "${__varname}=\"${__default}\""
    fi
}

# Validate ACCOUNT if provided on CLI against the account list
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
fi

# Interactive prompts for any options not supplied on CLI
if [[ -z "${TIME+x}" || -z "$TIME" ]]; then
    prompt_default TIME "Enter time (HH:MM:SS)" "$DEFAULT_TIME"
else
    echo "Using time: $TIME"
fi

if [[ -z "${NODES+x}" || -z "$NODES" ]]; then
    prompt_default NODES "Enter number of nodes" "$DEFAULT_NODES"
else
    echo "Using nodes: $NODES"
fi

if [[ -z "${CPUS+x}" || -z "$CPUS" ]]; then
    prompt_default CPUS "Enter number of CPUs per task" "$DEFAULT_CPUS"
else
    echo "Using cpus: $CPUS"
fi

if [[ -z "${MEM+x}" || -z "$MEM" ]]; then
    prompt_default MEM "Enter mem per CPU (e.g., 4096M)" "$DEFAULT_MEM"
else
    echo "Using mem: $MEM"
fi

# Account selection (interactive if not provided)
if [[ -n "$ACCOUNT" ]]; then
    echo "Using account: $ACCOUNT"
else
    if [[ $is_interactive -eq 1 ]]; then
        echo
        echo "Available $APERO_SERVER accounts:"
        echo
        for ((k=1;k<i;k++)); do
            echo "  $k) ${ACCOUNT_LIST[$k]}   (group: ${ACCOUNT_GROUP[$k]})"
        done
        if [[ -c /dev/tty ]]; then
            read -p "Select account by number: " CHOICE </dev/tty
        else
            read -p "Select account by number: " CHOICE
        fi
        ACCOUNT="${ACCOUNT_LIST[$CHOICE]}"
        if [[ -z "$ACCOUNT" ]]; then
            echo "Invalid selection."
            exit 1
        fi
    else
        # non-interactive: default to first account
        ACCOUNT="${ACCOUNT_LIST[1]}"
        echo "Non-interactive: defaulting account to: $ACCOUNT"
    fi
fi

# X11 handling: if WANT_X11 was set on CLI honor it; otherwise ask interactively
if [[ "$WANT_X11" == "Y" ]]; then
    X11_FLAG="--x11"
else
    if [[ $is_interactive -eq 1 ]]; then
        echo
        if [[ -c /dev/tty ]]; then
            read -p "Request an interactive X11 session? [Y/N]: " WANT_X11 </dev/tty
        else
            read -p "Request an interactive X11 session? [Y/N]: " WANT_X11
        fi
        WANT_X11=${WANT_X11:-N}
        if [[ "$WANT_X11" =~ ^[Yy]$ ]]; then
            X11_FLAG="--x11"
        else
            X11_FLAG=""
        fi
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

# Final action: if PROMPT=1 execute (we already asked); if PROMPT=0 do not run
# (interactive mode already collected options and displayed the command).
if [[ $PROMPT -eq 1 ]]; then
    if [[ $RUN_CONFIRMED -eq 1 ]]; then
        eval "$COMMAND"
    else
        # fallback - should not be reached, but ask just in case
        if [[ -c /dev/tty ]]; then
            read -p "Run salloc [Y]es or [N]o: " CONFIRM </dev/tty
        else
            read -p "Run salloc [Y]es or [N]o: " CONFIRM
        fi
        CONFIRM=${CONFIRM:-N}
        if [[ ! "$CONFIRM" =~ ^[Yy] ]]; then
            echo "Aborted by user."
            exit 1
        fi
        eval "$COMMAND"
    fi
else
    echo "No --prompt flag provided; command shown but not executed."
    echo "If you want to confirm and run the command add --prompt to the command line."
    exit 0
fi
