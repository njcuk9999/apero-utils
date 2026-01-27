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
INSTRUMENT_FILE="$APERO_BIN_PATH/apero_instruments.conf"
# Set the profile file
PROFILE_FILE="$APERO_BIN_PATH/apero_profiles.conf"
# User bashrc (used to detect whether apero_install added the install snippet)
USER_BASHRC="$HOME/.bashrc"

# -----------------------------------------------------------------------------
# Helper: get install_script for an instrument from INSTRUMENT_FILE
# returns relative path (as defined in conf) or empty string if not found
# -----------------------------------------------------------------------------
get_install_script_for_instrument() {
    local instr="$1"
    local current_section=""
    local line
    local install_script=""

    # guard: file exists
    if [[ ! -f "$INSTRUMENT_FILE" ]]; then
        echo ""
        return
    fi

    while IFS= read -r line || [[ -n "$line" ]]; do
        # trim whitespace (use bash parameter expansion)
        line="${line#${line%%[![:space:]]*}}"
        line="${line%${line##*[![:space:]]}}"
        # skip comments/empty
        [[ -z "$line" || "$line" =~ ^# ]] && continue

        if [[ "$line" =~ ^\[(.+)\]$ ]]; then
            current_section="${BASH_REMATCH[1]}"
            continue
        fi

        if [[ "$current_section" == "$instr" ]]; then
            if [[ "$line" =~ ^install_script[[:space:]]*=[[:space:]]*(.+)$ ]]; then
                install_script="${BASH_REMATCH[1]}"
                # return immediately
                echo "$install_script"
                return
            fi
        fi
    done < "$INSTRUMENT_FILE"

    echo ""
}

# -----------------------------------------------------------------------------
# Helper: is_instrument_installed
# - checks whether the install snippet (full path to install script) is present in ~/.bashrc
# - returns 0 if installed, 1 otherwise
# -----------------------------------------------------------------------------
is_instrument_installed() {
    local instr="$1"
    local rel_script
    rel_script=$(get_install_script_for_instrument "$instr")
    if [[ -z "$rel_script" ]]; then
        return 1
    fi
    local fullpath="$APERO_BIN_PATH/$rel_script"
    if [[ -f "$USER_BASHRC" ]] && grep -Fq "$fullpath" "$USER_BASHRC"; then
        return 0
    fi
    return 1
}

# -----------------------------------------------------------------------------
# Helper: list installed instruments (from INSTRUMENT_FILE that are present in ~/.bashrc)
# -----------------------------------------------------------------------------
list_installed_instruments() {
    if [[ ! -f "$INSTRUMENT_FILE" ]]; then
        echo "  (Instrument file not found: $INSTRUMENT_FILE)"
        return
    fi
    # read section names using awk to avoid sed issues
    awk '/^\[.*\]/{s=$0; gsub(/^\[|\]$/,"",s); print s}' "$INSTRUMENT_FILE" | while IFS= read -r instr; do
        if is_instrument_installed "$instr"; then
            echo "  - $instr"
        fi
    done
}

show_help() {
    echo "Usage: source activate.sh <instrument> <profile>"
    echo ""
    echo "This script must be sourced"
    echo ""
    echo "This script runs the commands defined for INSTRUMENT.PROFILE in apero_profiles.conf"
    echo ""
    echo "Available instruments (installed for this user):"
    if [[ -f "$INSTRUMENT_FILE" ]]; then
        list_installed_instruments
    else
        echo "  (Instrument file not found: $INSTRUMENT_FILE)"
    fi
    echo ""

    if [[ -n "$1" ]]; then
        local instr="$1"
        if [[ -f "$PROFILE_FILE" ]]; then
            # Use awk to safely extract profile names for this instrument
            local profiles
            profiles=$(awk -v inst="$instr" 'BEGIN{prefix="[" inst "."} index($0,prefix)==1 {s=substr($0,length(prefix)+1); gsub(/\].*$/,"",s); print s}' "$PROFILE_FILE")
            if [[ -n "$profiles" ]]; then
                echo "Available profiles for '$instr':"
                while IFS= read -r _p; do
                    printf '  - %s\n' "$_p"
                done <<< "$profiles"
            else
                echo "No profiles found for instrument '$instr'."
            fi
        else
            echo "Profile file not found: $PROFILE_FILE"
        fi
    fi

    echo ""
}

# -----------------------------------------------------------------------------
#  Detect if not sourced - show help and exit
# -----------------------------------------------------------------------------
(return 0 2>/dev/null)
if [ $? -ne 0 ]; then
    # If the script is executed rather than sourced, print a short, clear message
    # instead of the full help menu which can be noisy in automated contexts.
    echo "*************************"
    echo "ERROR: This script must be sourced, not executed."
    echo "Run: source $0 <instrument> <profile>"
    echo "*************************"
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
    echo "*************************"
    echo "ERROR: APERO_SERVER is not set."
    echo "Please run: $APERO_BIN_PATH/apero_install.sh"
    echo "*************************"
    return 1
fi

# -----------------------------------------------------------------------------
# Validate APERO_BIN_PATH and conf file
# -----------------------------------------------------------------------------
if [[ ! -f "$APERO_USERS_CONF" ]]; then
    echo "*************************"
    echo "ERROR: Cannot find apero_users.conf at:"
    echo "  $APERO_USERS_CONF"
    echo "Please contact the APERO administrators to be added."
    echo "*************************"
    return 1
fi

# -----------------------------------------------------------------------------
# Build lookup key [server.username]
# -----------------------------------------------------------------------------
LOOKUP="[$APERO_SERVER.$USER]"

# Find the line number where the header appears using awk (match at line start)
LINE=$(awk -v key="$LOOKUP" 'index($0,key)==1 {print NR; exit}' "$APERO_USERS_CONF")

# Check if header exists in the file
if [[ -z "$LINE" ]]; then
    echo "*************************"
    echo "ERROR: User entry '$LOOKUP' not found in apero_users.conf"
    echo "Please contact the APERO administrators to be added."
    echo "*************************"
    return 1
fi

# -----------------------------------------------------------------------------
# Extract name and email (lines after the header)
# -----------------------------------------------------------------------------
# name is next line, email the line after that
NAME=$(sed -n "$((LINE+1))p" "$APERO_USERS_CONF")
# trim name
NAME="${NAME#${NAME%%[![:space:]]*}}"
NAME="${NAME%${NAME##*[![:space:]]}}"
EMAIL=$(sed -n "$((LINE+2))p" "$APERO_USERS_CONF")
EMAIL="${EMAIL#${EMAIL%%[![:space:]]*}}"
EMAIL="${EMAIL%${EMAIL##*[![:space:]]}}"
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

# Detect --batch in the arguments passed when sourcing
BATCH=1
for _arg in "$@"; do
    if [[ "$_arg" == "--batch" ]]; then
        BATCH=0
        break
    fi
done

if [ -z "$INSTRUMENT" ]; then
    # Print a short message instead of the full help menu to avoid noisy output
    echo "*************************"
    echo "ERROR: No instrument selected."
    echo "*************************"
    echo ""
    echo "Available instruments (installed for this user):"
    if [[ -f "$INSTRUMENT_FILE" ]]; then
        list_installed_instruments
    else
        echo "  (Instrument file not found: $INSTRUMENT_FILE)"
    fi

    echo ""
    echo "Use: -h or --help for more details."
    echo ""
    return 1
fi


# -----------------------------------------------------------------------------
#  Validate instrument
# -----------------------------------------------------------------------------
# 1. Check if instrument actually exists in system instrument file
if ! grep -q "^\[$INSTRUMENT\]" "$INSTRUMENT_FILE"; then
    echo "*************************"
    echo "ERROR: Unknown instrument '$INSTRUMENT'"
    echo "*************************"
    echo ""
    echo "Available instruments (installed for this user):"
    if [[ -f "$INSTRUMENT_FILE" ]]; then
        list_installed_instruments
    else
        echo "  (Instrument file not found: $INSTRUMENT_FILE)"
    fi

    echo ""
    echo "Use: -h or --help for more details."
    echo ""
    return 1
fi

# 2. Require instrument to be installed via apero_install (present in ~/.bashrc)
if ! is_instrument_installed "$INSTRUMENT"; then
    echo "*************************"
    echo "ERROR: Instrument '$INSTRUMENT' exists but has not been installed for this user."
    echo "Please run: $APERO_BIN_PATH/apero_install.sh $INSTRUMENT"
    echo "*************************"
    echo ""
    echo "Available instruments (installed for this user):"
    if [[ -f "$INSTRUMENT_FILE" ]]; then
        list_installed_instruments
    else
        echo "  (Instrument file not found: $INSTRUMENT_FILE)"
    fi
    return 1
fi

# 3. Check if user is authorized to use the selected instrument
if [[ ! " $USER_INSTR_LIST " =~ " $INSTRUMENT " ]]; then
    echo "*************************"
    echo "ERROR: Instrument '$INSTRUMENT' exists but you are not authorized to use it."
    echo "Authorized instruments for $NAME: $USER_INSTR_LIST"
    echo "Please contact the APERO administrators to be added."
    echo "*************************"
    return 1
fi

# -----------------------------------------------------------------------------
#  Validate profile belongs to instrument
#    profile section looks like:  [nirps.profile1.v07]
# -----------------------------------------------------------------------------
FULL_SECTION="[$INSTRUMENT.$PROFILE]"

# Get list of profiles for this instrument
PROFILE_LIST=$(awk -v inst="$INSTRUMENT" 'BEGIN{prefix="[" inst "."} index($0,prefix)==1 {s=substr($0,length(prefix)+1); gsub(/\].*$/,"",s); print s}' "$PROFILE_FILE")

# -----------------------------------------------------------------------------
# CASE 1 — No profile provided
# -----------------------------------------------------------------------------
if [ -z "$PROFILE" ]; then
    echo "*************************"
    echo "ERROR: No profile selected for instrument '$INSTRUMENT'."
    echo "*************************"
    echo ""
    echo "Available profiles for '$INSTRUMENT':"
    if [[ -f "$PROFILE_FILE" ]]; then
        profiles=$(awk -v inst="$INSTRUMENT" 'BEGIN{prefix="[" inst "."} index($0,prefix)==1 {s=substr($0,length(prefix)+1); gsub(/\].*$/,"",s); print s}' "$PROFILE_FILE")
        if [[ -n "$profiles" ]]; then
            while IFS= read -r _p; do
                printf '  - %s\n' "$_p"
            done <<< "$profiles"
        else
            echo "  (No profiles found)"
        fi
    else
        echo "  (Profile file not found: $PROFILE_FILE)"
    fi

    echo ""
    echo "Use: -h or --help for more details."
    echo ""
    return 1
fi


# -----------------------------------------------------------------------------
# CASE 2 — Profile does not exist
# -----------------------------------------------------------------------------
if ! grep -q "^\[$INSTRUMENT\.$PROFILE\]" "$PROFILE_FILE"; then
    echo "*************************"
    echo "ERROR: Profile '$PROFILE' not found for instrument '$INSTRUMENT'."
    echo "*************************"
    echo ""
    echo "Available profiles for '$INSTRUMENT':"
    if [[ -f "$PROFILE_FILE" ]]; then
        profiles=$(awk -v inst="$INSTRUMENT" 'BEGIN{prefix="[" inst "."} index($0,prefix)==1 {s=substr($0,length(prefix)+1); gsub(/\].*$/,"",s); print s}' "$PROFILE_FILE")
        if [[ -n "$profiles" ]]; then
            while IFS= read -r _p; do
                printf '  - %s\n' "$_p"
            done <<< "$profiles"
        else
            echo "  (No profiles found)"
         fi
    else
        echo "  (Profile file not found: $PROFILE_FILE)"
    fi

    echo ""
    echo "Use: -h or --help for more details."
    echo ""
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
#  Add insturment-specific aliases
# -----------------------------------------------------------------------------
export APERO_INSTRUMENT="$INSTRUMENT"
# set the apero script path
APERO_SCRIPT_PATH="$APERO_PATH/${APERO_INSTRUMENT}_bin/scripts"
# software aliases
alias dfits="$APERO_SCRIPT_PATH/fitsio/dfits"
alias fitsort="$APERO_SCRIPT_PATH/fitsio/fitsort"
alias glow="$APERO_SCRIPT_PATH/glow/glow"

# apero tools
alias apero-trigger="cd $APERO_PATH/${APERO_INSTRUMENT}_bin/scripts/apero-utils/nirps/manual_trigger"
alias apero-checks="cd $APERO_PATH/${APERO_INSTRUMENT}_bin/scripts/apero-utils/nirps/apero_check"

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
echo "  gobin          : cd to the $INSTRUMENT bin directory"
echo "  godata         : cd to the $INSTRUMENT data directory"
echo "  dfits          : run dfits for $INSTRUMENT"
echo "  fitsort        : run fitsort for $INSTRUMENT"
echo "  glow           : run glow markdown viewer"
echo "  apero-trigger  : cd to manual trigger scripts for $INSTRUMENT"
echo "  apero-checks   : cd to APERO checks for $INSTRUMENT"
echo "  apero-activate : source the APERO profile activation script"
echo "  apero-salloc   : run APERO salloc launcher"
echo "  apero-find     : run APERO file finder tool"
echo ""
while IFS= read -r cmd; do
    eval "$cmd"
done < <(get_profile_commands "$INSTRUMENT.$PROFILE")

# If --batch was provided when sourcing, run the batch salloc command with
# the requested fixed options. This runs after profile commands have been
# executed.
if [[ "$BATCH" -eq 1 ]]; then
    SALLOC_SCRIPT="$APERO_BIN_PATH/apero_salloc.sh"
    echo "================================================="
    echo "Batch mode requested: launching salloc with preset options"
    echo "  Command: $SALLOC_SCRIPT --prompt --time=8 --cpus=2 --nodes=1 --mem=4096 --account=rrg-rdoyon"
    echo "================================================="
    if [[ -f "$SALLOC_SCRIPT" && -x "$SALLOC_SCRIPT" ]] || [[ -f "$SALLOC_SCRIPT" ]]; then
        # run the script (it will prompt/confirm as implemented in the salloc script)
        "$SALLOC_SCRIPT" --prompt --time=8 --cpus=2 --nodes=1 --mem=4096 --account=rrg-rdoyon
    else
        echo "*************************"
        echo "ERROR: salloc script not found at: $SALLOC_SCRIPT"
        echo "Skipping batch salloc launch."
        echo "*************************"
    fi
fi

