show_help() {
    cat <<EOF
Usage: apero_find.sh [START_PATH] [PATTERN] [COMMAND]

Find files matching a pattern and optionally run a command on each.

Arguments:
  START_PATH   Directory to start searching from.
  PATTERN      Filename glob pattern (default: "*").
  COMMAND      Command to run on each file (default: list only).

Examples:
  apero_find.sh /data
  apero_find.sh /data "*.fits"
  apero_find.sh /data "*.fits" "ls -l"
  apero_find.sh                   # fully interactive

Options:
  -h, --help    Show this help message and exit.
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
# 1. Get START_PATH
# -----------------------------------------------------------------------------
if [[ -n "$1" ]]; then
    START_PATH="$1"
else
    echo -n "Enter the starting path: "
    read START_PATH
fi

# Validate
if [[ ! -d "$START_PATH" ]]; then
    echo "Error: '$START_PATH' is not a directory."
    exit 1
fi

cd "$START_PATH" || exit 1

# -----------------------------------------------------------------------------
# 2. Get PATTERN
# -----------------------------------------------------------------------------
if [[ -n "$2" ]]; then
    PATTERN="$2"
else
    echo -n "Enter filename pattern (blank = *): "
    read PATTERN
fi

# Default to "*"
if [[ -z "$PATTERN" ]]; then
    PATTERN="*"
fi

echo "Using pattern: $PATTERN"

# -----------------------------------------------------------------------------
# 3. Get COMMAND
# -----------------------------------------------------------------------------
if [[ -n "$3" ]]; then
    USER_CMD="$3"
else
    echo -n "Enter command to run on each file (blank to list only): "
    read USER_CMD
fi

echo
echo "Searching in: $(pwd)"
echo

# -----------------------------------------------------------------------------
# 4. Execute the find
# -----------------------------------------------------------------------------
if [[ -z "$USER_CMD" ]]; then
    # List files only
    find "$(pwd)" -type f -name "$PATTERN"
else
    # Run a command on each file
    find "$(pwd)" -type f -name "$PATTERN" -print0 | \
    while IFS= read -r -d '' FILE; do
        eval "$USER_CMD \"$FILE\""
    done
fi