#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# CODE NAME HERE

# CODE DESCRIPTION HERE

Created on 2026-01-29 at 09:05

@author: cook
"""
from typing import List
import os
import argparse


# =============================================================================
# Define variables
# =============================================================================
# valid extensions
VALID_EXTENSIONS = ['.fits']


# =============================================================================
# Define functions
# =============================================================================
def get_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Index and compare raw file lists')
    parser.add_argument('--index', type=str, default=None,
                        help='Path to raw directory to index')
    parser.add_argument('--savefile', type=str, default='raw_files.csv',
                        help='Path to save the index file')
    parser.add_argument('--list1', type=str, default=None,
                        help='First list file to compare')
    parser.add_argument('--list2', type=str, default=None,
                        help='Second list file to compare')
    return parser.parse_args()


def get_files(raw_path: str) -> List[str]:
    """
    Recursively find all files with valid extensions in raw_path

    :param raw_path: Path to raw data directory
    :return: List of relative file paths
    """
    files = []

    # Expand user path if needed
    raw_path = os.path.expanduser(raw_path)

    # Check if path exists
    if not os.path.exists(raw_path):
        print(f"Error: Path {raw_path} does not exist")
        return files

    # Walk through directory recursively
    print(f"Scanning directory: {raw_path}")
    for root, dirs, filenames in os.walk(raw_path):
        for filename in filenames:
            # Check if file has valid extension
            _, ext = os.path.splitext(filename)
            if ext.lower() in VALID_EXTENSIONS:
                # Get full path
                full_path = os.path.join(root, filename)
                # Get relative path from raw_path
                rel_path = os.path.relpath(full_path, raw_path)
                files.append(rel_path)

    print(f"Found {len(files)} files with valid extensions")
    return files


def save_index(save_file: str, file_list: List[str]):
    """
    Save file list to CSV file

    :param save_file: Path to save file
    :param file_list: List of file paths to save
    """
    # Expand user path if needed
    save_file = os.path.expanduser(save_file)

    # Create directory if it doesn't exist
    save_dir = os.path.dirname(save_file)
    if save_dir and not os.path.exists(save_dir):
        os.makedirs(save_dir)

    # Save to CSV
    with open(save_file, 'w') as f:
        for file_path in file_list:
            f.write(f"{file_path}\n")

    print(f"Saved {len(file_list)} files to {save_file}")


def compare_files(list1_path: str, list2_path: str) -> List[str]:
    """
    Compare two file lists and return unmatched files

    :param list1_path: Path to first list file
    :param list2_path: Path to second list file
    :return: List of unmatched file paths
    """
    # Expand user paths if needed
    list1_path = os.path.expanduser(list1_path)
    list2_path = os.path.expanduser(list2_path)

    # Check if files exist
    if not os.path.exists(list1_path):
        print(f"Error: File {list1_path} does not exist")
        return []
    if not os.path.exists(list2_path):
        print(f"Error: File {list2_path} does not exist")
        return []

    # Read file lists
    with open(list1_path, 'r') as f:
        list1 = set(line.strip() for line in f if line.strip())

    with open(list2_path, 'r') as f:
        list2 = set(line.strip() for line in f if line.strip())

    # Find files that are in list1 but not in list2, and vice versa
    only_in_list1 = list1 - list2
    only_in_list2 = list2 - list1
    matched = list1 & list2

    # Print statistics
    print(f"\n{'='*60}")
    print(f"Comparison Results:")
    print(f"{'='*60}")
    print(f"Files in {os.path.basename(list1_path)}: {len(list1)}")
    print(f"Files in {os.path.basename(list2_path)}: {len(list2)}")
    print(f"Files matched: {len(matched)}")
    print(f"Files only in {os.path.basename(list1_path)}: {len(only_in_list1)}")
    print(f"Files only in {os.path.basename(list2_path)}: {len(only_in_list2)}")
    print(f"Total unmatched files: {len(only_in_list1) + len(only_in_list2)}")
    print(f"{'='*60}\n")

    # Return all unmatched files (sorted)
    unmatched = sorted(list(only_in_list1) + list(only_in_list2))
    return unmatched


def main():
    # ----------------------------------------------------------------------
    # get arguments:
    #   --index={raw_path} to index a raw directory
    #   --savefile={file.csv} to save the index to
    #   --list1={file.csv} to compare against list2
    #   --list2={file.csv} to compare against list1
    args = get_args()

    if args.index not in ['', 'None', None]:
        raw_files = get_files(raw_path=args.index)

        save_index(args.savefile, raw_files)
        return

    # otherwise we should compare two lists
    cond1 = args.list1 not in ['', 'None', None]
    cond2 = args.list2 not in ['', 'None', None]
    if cond1 and cond2:
        # compare two lists
        unmatched_files = compare_files(list1_path=args.list1,
                                        list2_path=args.list2)
        if len(unmatched_files) > 0:
            # generate a compfile name
            part1 = os.path.basename(args.list1).replace('.csv', '_vs_')
            part2 = os.path.basename(args.list2).replace('.csv',
                                                         '_unmatched.csv')
            compfile = str(part1) + str(part2)

            save_index(compfile, unmatched_files)
        else:
            print("All files matched! No unmatched files to save.")
        return

    # If neither indexing nor comparing, show help
    print("Error: Please specify either --index or both --list1 and --list2")
    print("Use --help for more information")


# =============================================================================
# Start of code
# =============================================================================
# Main code here
if __name__ == "__main__":
    # run main function
    main()

# =============================================================================
# End of code
# =============================================================================
