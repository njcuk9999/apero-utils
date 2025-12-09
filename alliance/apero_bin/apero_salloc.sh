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

# Time
read -p "Enter time (HH:MM:SS) [default 4:00:00]: " TIME
TIME=${TIME:-4:00:00}

# Nodes
read -p "Enter number of nodes [default 1]: " NODES
NODES=${NODES:-1}

# CPUs
read -p "Enter number of CPUs per task [default 20]: " CPUS
CPUS=${CPUS:-20}

# Memory
read -p "Enter mem per CPU (e.g., 4096M) [default 4096M]: " MEM
MEM=${MEM:-4096M}

echo
echo "Available APERO accounts:"
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

echo
echo "Running:"
echo "salloc --time=$TIME --cpus-per-task=$CPUS --nodes=$NODES --mem-per-cpu=$MEM --account=$ACCOUNT"
echo