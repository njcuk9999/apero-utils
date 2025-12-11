# These need changing on a different server or project
export APERO_PROJECT_ID="6102120"
export APERO_SERVER="alliance"
TEST_BIN_PATH="/project/$APERO_PROJECT_ID/apero/"

# -----------------------------------------------------------------------------
# Do not change under here
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Function: Show help message
# -----------------------------------------------------------------------------
show_help() {
    echo ""
    echo "Do not use - this script is just for sourcing elsewhere"
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
# Determine APERO_BIN_PATH
if [[ -n "$APERO_PROJECT_ID" && -d $TEST_BIN_PATH ]]; then
    APERO_PATH=$TEST_BIN_PATH
else
    # Get the directory of this script
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    APERO_PATH="$(cd "$(dirname "${SCRIPT_DIR}")" && pwd)"
fi
# set the bin path
APERO_BIN_PATH="$APERO_PATH/apero_bin"

export APERO_PATH=$PATH_PATH
export APERO_BIN_PATH=$APERO_BIN_PATH