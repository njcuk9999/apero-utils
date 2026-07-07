#!/bin/bash

# Define directories
dir1="/cosmos99/nirps/git-bin/apero-utils/nirps/manual_trigger/.apero/manual_trigger/"
dir2="/home/nirps-client/.apero/manual_trigger/"
backup_dir="/cosmos99/nirps/git-bin/apero-utils/nirps/manual_trigger/.apero/manual_trigger_backup_$(date +%Y%m%d_%H%M%S)"

# 1. Back up both directories
echo "Backing up directories..."
mkdir -p "$backup_dir"
cp -r "$dir1" "$backup_dir/dir1"
cp -r "$dir2" "$backup_dir/dir2"
echo "Backup created at: $backup_dir"

# 2. Merge logs
# We iterate through files in dir2 and merge them into dir1
for file in "$dir2"/*.log; do
    filename=$(basename "$file")
    target="$dir1/$filename"

    if [ -f "$target" ]; then
        echo "Merging $filename..."
        # Using awk to keep lines from dir1, then adding unique lines from dir2
        # !seen[$0]++ ensures we don't add duplicate lines
        awk '!seen[$0]++' "$target" "$file" > "$target.tmp"
        mv "$target.tmp" "$target"
    else
        echo "File $filename not found in $dir1, skipping."
    fi
done

echo "Merge complete."