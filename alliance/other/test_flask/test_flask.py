#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Flask File Explorer Test App

A simple Flask web application that provides an in-browser file explorer
for a specified directory and its subdirectories.

Created on 2026-01-20 at 13:10

@author: cook
"""

import os
from flask import Flask, render_template, abort, request
import argparse

# =============================================================================
# Define variables
# =============================================================================
PATH_TO_FILES = '~/ari/data/'
DEFAULT_PORT = 5000
# -----------------------------------------------------------------------------

# =============================================================================
# HTML Template
# =============================================================================
HTML_TEMPLATE = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>ARI Flask Test</title>
    <style>
        body {
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            max-width: 1200px;
            margin: 0 auto;
            padding: 20px;
            background-color: #f5f5f5;
        }
        h1 {
            color: #2c3e50;
            border-bottom: 3px solid #3498db;
            padding-bottom: 10px;
        }
        .info-box {
            background-color: white;
            border-radius: 8px;
            padding: 20px;
            margin: 20px 0;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }
        .info-item {
            margin: 10px 0;
            padding: 8px;
            background-color: #ecf0f1;
            border-left: 4px solid #3498db;
            border-radius: 4px;
        }
        .info-label {
            font-weight: bold;
            color: #2c3e50;
        }
        .explorer {
            background-color: white;
            border-radius: 8px;
            padding: 20px;
            margin: 20px 0;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }
        .breadcrumb {
            background-color: #ecf0f1;
            padding: 10px;
            border-radius: 4px;
            margin-bottom: 15px;
            font-family: monospace;
        }
        .file-list {
            list-style: none;
            padding: 0;
        }
        .file-item, .dir-item {
            padding: 10px;
            margin: 5px 0;
            border-radius: 4px;
            transition: background-color 0.2s;
        }
        .file-item:hover, .dir-item:hover {
            background-color: #ecf0f1;
        }
        .dir-item {
            background-color: #e8f4f8;
        }
        .dir-item a {
            color: #2980b9;
            text-decoration: none;
            font-weight: bold;
        }
        .dir-item a:hover {
            text-decoration: underline;
        }
        .file-item {
            background-color: #f9f9f9;
            color: #555;
        }
        .icon {
            margin-right: 8px;
            font-size: 1.2em;
        }
        .parent-link {
            background-color: #fff3cd;
            border-left: 4px solid #ffc107;
        }
        .stats {
            display: flex;
            gap: 20px;
            margin-top: 10px;
        }
        .stat-box {
            flex: 1;
            padding: 15px;
            background-color: #3498db;
            color: white;
            border-radius: 4px;
            text-align: center;
        }
        .stat-value {
            font-size: 2em;
            font-weight: bold;
        }
        .stat-label {
            font-size: 0.9em;
            opacity: 0.9;
        }
    </style>
</head>
<body>
    <h1>🗂️ ARI Flask Test</h1>
    
    <div class="info-box">
        <h2>Server Information</h2>
        <div class="info-item">
            <span class="info-label">Port Number:</span> {{ port }}
        </div>
        <div class="info-item">
            <span class="info-label">Root Path:</span> <code>{{ root_path }}</code>
        </div>
        <div class="info-item">
            <span class="info-label">Current Path:</span> <code>{{ current_path }}</code>
        </div>
        <div class="stats">
            <div class="stat-box">
                <div class="stat-value">{{ dir_count }}</div>
                <div class="stat-label">Directories</div>
            </div>
            <div class="stat-box" style="background-color: #2ecc71;">
                <div class="stat-value">{{ file_count }}</div>
                <div class="stat-label">Files</div>
            </div>
            <div class="stat-box" style="background-color: #9b59b6;">
                <div class="stat-value">{{ total_count }}</div>
                <div class="stat-label">Total Items</div>
            </div>
        </div>
    </div>
    
    <div class="explorer">
        <h2>📁 File Explorer</h2>
        <div class="breadcrumb">
            {{ current_path }}
        </div>
        
        <ul class="file-list">
            {% if show_parent %}
            <li class="dir-item parent-link">
                <span class="icon">⬆️</span>
                <a href="?path={{ parent_path }}">.. (Parent Directory)</a>
            </li>
            {% endif %}
            
            {% for item in directories %}
            <li class="dir-item">
                <span class="icon">📁</span>
                <a href="?path={{ item.path }}">{{ item.name }}</a>
            </li>
            {% endfor %}
            
            {% for item in files %}
            <li class="file-item">
                <span class="icon">📄</span>
                {{ item.name }} <span style="color: #999; font-size: 0.9em;">({{ item.size }})</span>
            </li>
            {% endfor %}
            
            {% if not directories and not files %}
            <li style="padding: 20px; text-align: center; color: #999;">
                <em>This directory is empty</em>
            </li>
            {% endif %}
        </ul>
    </div>
</body>
</html>
"""

# =============================================================================
# Define functions
# =============================================================================
def format_size(size_bytes):
    """Format file size in human-readable format"""
    for unit in ['B', 'KB', 'MB', 'GB', 'TB']:
        if size_bytes < 1024.0:
            return f"{size_bytes:.1f} {unit}"
        size_bytes /= 1024.0
    return f"{size_bytes:.1f} PB"


def get_directory_contents(path):
    """Get directories and files in the given path"""
    directories = []
    files = []

    try:
        entries = sorted(os.listdir(path))

        for entry in entries:
            full_path = os.path.join(path, entry)

            try:
                if os.path.isdir(full_path):
                    directories.append({
                        'name': entry,
                        'path': full_path
                    })
                elif os.path.isfile(full_path):
                    size = os.path.getsize(full_path)
                    files.append({
                        'name': entry,
                        'size': format_size(size)
                    })
            except (PermissionError, OSError):
                # Skip files/dirs we can't access
                continue

    except (PermissionError, OSError) as e:
        print(f"Error accessing directory: {e}")

    return directories, files


def create_app(root_path, port):
    """Create and configure Flask app"""
    app = Flask(__name__)

    # Expand the root path to handle ~
    root_path = os.path.expanduser(root_path)
    root_path = os.path.abspath(root_path)

    # Create directory if it doesn't exist
    if not os.path.exists(root_path):
        print(f"Warning: Directory '{root_path}' does not exist. Creating it...")
        os.makedirs(root_path, exist_ok=True)

    @app.route('/')
    def index():

        # Get the requested path from query parameter
        requested_path = request.args.get('path', root_path)
        requested_path = os.path.abspath(requested_path)

        # Security check: ensure we're not accessing outside root_path
        if not requested_path.startswith(root_path):
            abort(403, "Access denied: Cannot access paths outside root directory")

        if not os.path.exists(requested_path):
            abort(404, "Path not found")

        if not os.path.isdir(requested_path):
            abort(400, "Path is not a directory")

        # Get directory contents
        directories, files = get_directory_contents(requested_path)

        # Determine if we should show parent directory link
        show_parent = requested_path != root_path
        parent_path = os.path.dirname(requested_path) if show_parent else None

        # Count totals
        dir_count = len(directories)
        file_count = len(files)
        total_count = dir_count + file_count

        return render_template(
            'index.html',
            port=port,
            root_path=root_path,
            current_path=requested_path,
            directories=directories,
            files=files,
            dir_count=dir_count,
            file_count=file_count,
            total_count=total_count,
            show_parent=show_parent,
            parent_path=parent_path
        )

    return app


# =============================================================================
# Start of code
# =============================================================================
if __name__ == "__main__":
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='ARI Flask File Explorer')
    parser.add_argument('--path', type=str, default=PATH_TO_FILES,
                        help=f'Path to explore (default: {PATH_TO_FILES})')
    parser.add_argument('--port', type=int, default=DEFAULT_PORT,
                        help=f'Port number (default: {DEFAULT_PORT})')
    parser.add_argument('--host', type=str, default='127.0.0.1',
                        help='Host address (default: 127.0.0.1)')

    args = parser.parse_args()

    # Create and run the app
    app = create_app(args.path, args.port)

    print("=" * 60)
    print("ARI Flask File Explorer")
    print("=" * 60)
    print(f"Root Path: {os.path.expanduser(args.path)}")
    print(f"Host: {args.host}")
    print(f"Port: {args.port}")
    print(f"URL: http://{args.host}:{args.port}")
    print("=" * 60)
    print("Press Ctrl+C to stop the server")
    print("=" * 60)

    app.run(host=args.host, port=args.port, debug=True)

# =============================================================================
# End of code
# =============================================================================
