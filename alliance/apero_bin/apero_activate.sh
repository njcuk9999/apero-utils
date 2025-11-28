#!/bin/bash

#  This script runs the commands defined for INSTRUMENT.PROFILE in apero_profiles.conf

# -----------------------------
#  Detect if sourced
# -----------------------------
(return 0 2>/dev/null)
if [ $? -ne 0 ]; then
    echo "ERROR: This script must be sourced, not executed."
    echo "Use: source activate.sh <instrument> <profile>"
    exit 1
fi

# -----------------------------
#  Basic setup & arguments
# -----------------------------
INSTRUMENT="$1"
PROFILE="$2"

if [ -z "$INSTRUMENT" ]; then
    echo "Usage: source activate.sh <instrument> <profile>"
    echo ""
    echo "Available instruments:"
    cut -d= -f1 "$APERO_BIN_PATH/apero_instruments.ini"
    return 1
fi

APERO_BIN_PATH="/project/$APERO_PROJECT_ID/apero/apero_bin"
INSTRUMENT_FILE="$APERO_BIN_PATH/apero_instruments.ini"
PROFILE_FILE="$APERO_BIN_PATH/apero_profiles.conf"

# -----------------------------
#  Validate instrument
# -----------------------------
if ! grep -q "^$INSTRUMENT=" "$INSTRUMENT_FILE"; then
    echo "ERROR: Unknown instrument '$INSTRUMENT'"
    echo "Available instruments:"
    cut -d= -f1 "$INSTRUMENT_FILE"
    return 1
fi


# -----------------------------
#  Validate profile belongs to instrument
#    profile section looks like:  [nirps.profile1.v07]
# -----------------------------
FULL_SECTION="[$INSTRUMENT.$PROFILE]"

# Get list of profiles for this instrument
PROFILE_LIST=$(grep "^\[$INSTRUMENT\." "$PROFILE_FILE" \
                | sed "s/^\[$INSTRUMENT\.//; s/\].*$//")

# -------------------------------
# CASE 1 — No profile provided
# -------------------------------
if [ -z "$PROFILE" ]; then
    echo "No profile selected for instrument '$INSTRUMENT'."
    echo "Available profiles:"
    echo "    $PROFILE_LIST"
    return 1
fi

# -------------------------------
# CASE 2 — Profile does not exist
# -------------------------------
if ! grep -q "^\[$INSTRUMENT\.$PROFILE\]" "$PROFILE_FILE"; then
    echo "ERROR: Profile '$PROFILE' not found for instrument '$INSTRUMENT'."
    echo "Available profiles:"
    echo "    $PROFILE_LIST"
    return 1
fi


# -----------------------------
#  Function: Read commands in a section
# -----------------------------
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


# -----------------------------
#  Run commands in profile
# -----------------------------
echo "================================================="
echo "Welcome to APERO-$INSTRUMENT @ Alliance"
echo "================================================="
echo "Activating instrument: $INSTRUMENT"
echo "Using profile:        $PROFILE"
echo "================================================="
echo ""
echo ""

while IFS= read -r cmd; do
    eval "$cmd"
done < <(get_profile_commands "$INSTRUMENT.$PROFILE")

