#!/bin/bash

#  This script runs the commands defined for INSTRUMENT.PROFILE in apero_profiles.conf

# -----------------------------------------------------------------------------
# set up variables
# -----------------------------------------------------------------------------
# Get core source script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source $SCRIPT_DIR/apero_core.sh
# Set the user configuration file
APERO_USERS_CONF="$APERO_BIN_PATH/apero_users.conf"
# Set the instrument file
INSTRUMENT_FILE="$APERO_BIN_PATH/apero_instruments.ini"
# Set the porfile file
PROFILE_FILE="$APERO_BIN_PATH/apero_profiles.conf"

show_help() {
    echo "Usage: source activate.sh <instrument> <profile>"
    echo ""
    echo "This script must be sourced"
    echo ""
    echo "This script runs the commands defined for INSTRUMENT.PROFILE in apero_profiles.conf"
    echo ""
    echo "Available instruments:"
    if [[ -f "$INSTRUMENT_FILE" ]]; then
        cut -d= -f1 "$INSTRUMENT_FILE" | sed 's/^/  - /'
    else
        echo "  (Instrument file not found: $INSTRUMENT_FILE)"
    fi
    echo ""

    if [[ -n "$1" ]]; then
        local instr="$1"
        if [[ -f "$PROFILE_FILE" ]]; then
            local profiles=$(grep "^\[$instr\." "$PROFILE_FILE" | sed "s/^\[$instr\.//; s/\].*$//")
            if [[ -n "$profiles" ]]; then
                echo "Available profiles for '$instr':"
                echo "$profiles" | sed 's/^/  - /'
            else
                echo "No profiles found for instrument '$instr'."
            fi
        else
            echo "Profile file not found: $PROFILE_FILE"
        fi
    fi

    echo ""
    echo "Example:"
    echo "  source activate.sh nirps profile1"
    echo ""
}

# -----------------------------------------------------------------------------
#  Detect if not sourced - show help and exit
# -----------------------------------------------------------------------------
(return 0 2>/dev/null)
if [ $? -ne 0 ]; then
    show_help
    exit 1
fi

# -----------------------------------------------------------------------------
#  Check for help flag
# -----------------------------------------------------------------------------
if [[ "$1" == "-h" || "$1" == "--help" ]]; then
    show_help
    return 0
fi

# -----------------------------------------------------------------------------
# Check APERO_SERVER is set
# -----------------------------------------------------------------------------
if [[ -z "$APERO_SERVER" ]]; then
    echo "ERROR: APERO_SERVER is not set."
    echo "Please run: $APERO_BIN_PATH/apero_install.sh"
    return 1
fi

# -----------------------------------------------------------------------------
# Validate APERO_BIN_PATH and conf file
# -----------------------------------------------------------------------------
if [[ ! -f "$APERO_USERS_CONF" ]]; then
    echo "ERROR: Cannot find apero_users.conf at:"
    echo "  $APERO_USERS_CONF"
    echo "Please contact the APERO administrators to be added."
    return 1
fi

# -----------------------------------------------------------------------------
# Build lookup key [server.username]
# -----------------------------------------------------------------------------
LOOKUP="[$APERO_SERVER.$USER]"
ESCAPED_LOOKUP=$(printf '%s\n' "$LOOKUP" | sed 's/[][\.^$*+?{|}()]/\\&/g')

# Check if header exists in the file
if ! grep -q "^$ESCAPED_LOOKUP" "$APERO_USERS_CONF"; then
    echo "ERROR: User entry '$LOOKUP' not found in apero_users.conf"
    echo "Please contact the APERO administrators to be added."
    return 1
fi

# -----------------------------------------------------------------------------
# Extract name and email (lines after the header)
# -----------------------------------------------------------------------------
# Get the line number where the header appears
LINE=$(grep -n "^$ESCAPED_LOOKUP" "$APERO_USERS_CONF" | head -n 1 | cut -d: -f1)

# ensure LINE is numeric
if ! [[ "$LINE" =~ ^[0-9]+$ ]]; then
    echo "ERROR: Could not locate user header '$LOOKUP' in $APERO_USERS_CONF"
    echo "Please contact the APERO administrators to be added."
    return 1
fi

# name is next line, email the line after that
NAME=$(sed -n "$((LINE+1))p" "$APERO_USERS_CONF")
EMAIL=$(sed -n "$((LINE+2))p" "$APERO_USERS_CONF")
USER_INSTR=$(sed -n "$((LINE+3))p" "$APERO_USERS_CONF" | tr -d '[:space:]')
# Convert comma list → space list
USER_INSTR_LIST=$(echo "$USER_INSTR" | tr ',' ' ')

# -----------------------------------------------------------------------------
# Export environment variables for user
# -----------------------------------------------------------------------------
export APERO_USER="$USER"
export APERO_USER_NAME="$NAME"
export APERO_USER_EMAIL="$EMAIL"

# -----------------------------------------------------------------------------
#  Basic setup & arguments for apero-activate
# -----------------------------------------------------------------------------
INSTRUMENT="$1"
PROFILE="$2"

if [ -z "$INSTRUMENT" ]; then
    show_help
    return 1
fi


# -----------------------------------------------------------------------------
#  Validate instrument
# -----------------------------------------------------------------------------
# 1. Check if instrument actually exists in system instrument file
if ! grep -q "^$INSTRUMENT=" "$INSTRUMENT_FILE"; then
    echo "ERROR: Unknown instrument '$INSTRUMENT'"
    echo "Available instruments:"
    cut -d= -f1 "$INSTRUMENT_FILE"
    return 1
fi
# 2. Check if user is authorized to use the selected instrument
if [[ ! " $USER_INSTR_LIST " =~ " $INSTRUMENT " ]]; then
    echo "ERROR: Instrument '$INSTRUMENT' exists but you are not authorized to use it."
    echo "Authorized instruments for $NAME: $USER_INSTR_LIST"
    echo "Please contact the APERO administrators to be added."
    return 1
fi

# -----------------------------------------------------------------------------
#  Validate profile belongs to instrument
#    profile section looks like:  [nirps.profile1.v07]
# -----------------------------------------------------------------------------
FULL_SECTION="[$INSTRUMENT.$PROFILE]"

# Get list of profiles for this instrument
PROFILE_LIST=$(grep "^\[$INSTRUMENT\." "$PROFILE_FILE" \
                | sed "s/^\[$INSTRUMENT\.//; s/\].*$//")

# -----------------------------------------------------------------------------
# CASE 1 — No profile provided
# -----------------------------------------------------------------------------
if [ -z "$PROFILE" ]; then
    echo "No profile selected for instrument '$INSTRUMENT'."
    show_help "$INSTRUMENT"
    return 1
fi


# -----------------------------------------------------------------------------
# CASE 2 — Profile does not exist
# -----------------------------------------------------------------------------
if ! grep -q "^\[$INSTRUMENT\.$PROFILE\]" "$PROFILE_FILE"; then
    echo "ERROR: Profile '$PROFILE' not found for instrument '$INSTRUMENT'."
    show_help "$INSTRUMENT"
    return 1
fi

# -----------------------------------------------------------------------------
#  Function: Read commands in a section
# -----------------------------------------------------------------------------
get_profile_commands() {
    local section="[$1]"
    local in_section=0

    while IFS= read -r line; do

        # start of section
        if [[ "$line" == "$section" ]]; then
            in_section=1
            continue
        fi

        # entering new section → stop
        if [[ "$line" =~ ^\[.*\]$ ]]; then
            [ $in_section -eq 1 ] && break
        fi

        # return only non-empty, non-comment lines
        if [ $in_section -eq 1 ]; then
            [[ -z "$line" ]] && continue
            [[ "$line" =~ ^# ]] && continue
            echo "$line"
        fi

    done < "$PROFILE_FILE"
}


# -----------------------------------------------------------------------------
#  Run commands in profile
# -----------------------------------------------------------------------------
echo "================================================="
echo "Welcome to APERO-$INSTRUMENT @ Alliance"
echo "================================================="
echo "User = $APERO_USER_NAME [$APERO_USER]"
echo "Email = $APERO_USER_EMAIL"
echo "-------------------------------------------------"
echo "Activating instrument: $INSTRUMENT"
echo "Using profile:        $PROFILE"
echo "================================================="
echo ""
echo "Available aliases:"
echo "  goapero        : cd to the APERO project directory"
echo "  gobin          : cd to the $APERO_INSTRUMENT bin directory"
echo "  godata         : cd to the $APERO_INSTRUMENT data directory"
echo "  dfits          : run dfits for $INSTRUMENT"
echo "  fitsort        : run fitsort for $INSTRUMENT"
echo "  apero-trigger  : cd to manual trigger scripts for $INSTRUMENT"
echo "  apero-checks   : cd to APERO checks for $INSTRUMENT"
echo "  apero-activate : source the APERO profile activation script"
echo "  apero-salloc   : run APERO salloc launcher"
echo "  apero-find     : run APERO file finder tool"
echo ""

while IFS= read -r cmd; do
    eval "$cmd"
done < <(get_profile_commands "$INSTRUMENT.$PROFILE")

