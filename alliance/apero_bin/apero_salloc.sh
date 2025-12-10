#!/bin/bash

# define the apero bin path
APERO_BIN_PATH="/project/$APERO_PROJECT_ID/apero/apero_bin"

# define conf file
CONF_FILE="$APERO_BIN_PATH/apero_groups.conf"

# Test if conf file exists
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
# Ask about email
# -----------------------------------------------------------------------------
EMAIL_FLAG=""

# If APERO_USER_EMAIL exists, ask whether to use it
if [[ -n "$APERO_USER_EMAIL" ]]; then
    echo
    read -p "Use email notifications for $APERO_USER_EMAIL? [Y/n]: " USE_EMAIL
    USE_EMAIL=${USE_EMAIL:-Y}

    if [[ "$USE_EMAIL" =~ ^[Yy]$ ]]; then
        EMAIL_FLAG="--mail-type=ALL --mail-user=$APERO_USER_EMAIL"
    fi

else
    # Email variable NOT set → ask user to enter manually
    echo
    read -p "Enter email for notifications (leave blank for none): " ENTERED_EMAIL
    if [[ -n "$ENTERED_EMAIL" ]]; then
        EMAIL_FLAG="--mail-type=ALL --mail-user=$ENTERED_EMAIL"
    fi
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

COMMAND="salloc --time=$TIME --cpus-per-task=$CPUS --nodes=$NODES --mem-per-cpu=$MEM --account=$ACCOUNT $EMAIL_FLAG $X11_FLAG"

echo
echo "The following salloc command will be run:"
echo
echo " $COMMAND"
echo

read -p "Run this command? [Y/n]: " CONFIRM
CONFIRM=${CONFIRM:-Y}

if [[ ! "$CONFIRM" =~ ^[Yy]$ ]]; then
    echo "Aborted by user."
    exit 1
fi

eval "$COMMAND"