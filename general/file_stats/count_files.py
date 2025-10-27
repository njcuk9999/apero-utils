#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2025-10-27 at 13:32

@author: cook
"""
import os
from tqdm import tqdm

# =============================================================================
# Define variables
# =============================================================================


# =============================================================================
# Define functions
# =============================================================================
def count_files(base_dir: str):
    # First count how many files exist (for tqdm total)
    all_files = []

    pbar = tqdm(os.walk(base_dir, followlinks=False), desc='Indexing files...')

    for root, _, files in pbar:
        for f in files:
            all_files.append(os.path.join(root, f))
            # Dynamically update description
            pbar.set_description(f'Indexed {len(all_files)} files... [root={root}]')

    pbar.close()

    total_files = len(all_files)
    symlinks = 0
    non_symlinks = 0
    total_size = 0

    for path in tqdm(all_files, desc="Scanning files", unit="file"):
        try:
            if os.path.islink(path):
                symlinks += 1
            else:
                non_symlinks += 1
                total_size += os.path.getsize(path)
        except OSError:
            pass

    return {
        "total_files": total_files,
        "symlinks": symlinks,
        "non_symlinks": non_symlinks,
        "total_size_bytes": total_size
    }


# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":
    # ----------------------------------------------------------------------
    import sys
    path = sys.argv[1] if len(sys.argv) > 1 else "."
    stats = count_files(path)
    print(f"Directory: {os.path.abspath(path)}")
    print(f"Total files: {stats['total_files']}")
    print(f"Symlinks: {stats['symlinks']}")
    print(f"Non-symlinks: {stats['non_symlinks']}")
    print(f"Total size: {stats['total_size_bytes'] / 1_048_576:.2f} MB")

# =============================================================================
# End of code
# =============================================================================
