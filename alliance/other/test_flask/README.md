# ARI Flask Test - File Explorer

A simple Flask web application that provides an in-browser file explorer for browsing directories and subdirectories.

## Features

- 📁 Browse any directory on your filesystem
- 📊 Display server information (port, paths, file/directory counts)
- 🔍 Navigate through subdirectories with click-through interface
- 📏 Show file sizes in human-readable format
- 🔒 Security check to prevent access outside root directory
- 🎨 Clean, modern interface with visual distinction between files and folders

## Installation

### Option 1: Using Conda (Recommended)

```bash
# Create a new conda environment
conda create -n flask-test python=3.9

# Activate the environment
conda activate flask-test

# Install dependencies
pip install -r requirements.txt
```

### Option 2: Using venv (Python Virtual Environment)

```bash
# Create a virtual environment
python -m venv venv

# Activate the environment
# On Linux/Mac:
source venv/bin/activate
# On Windows:
# venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt
```

### Option 3: System-wide Installation (Not Recommended)

```bash
# Install directly (not recommended for production)
pip install -r requirements.txt
```

## Usage

### Basic Usage

Run the application with default settings (explores `~/ari/data/` on port 5000):

```bash
python test_flask.py
```

Then open your web browser and navigate to:
```
http://127.0.0.1:5000
```

### Custom Path

Specify a different directory to explore:

```bash
python test_flask.py --path /path/to/your/directory
```

### Custom Port

Use a different port:

```bash
python test_flask.py --port 8080
```

Then navigate to: `http://127.0.0.1:8080`

### Custom Host (Allow External Access)

Allow access from other machines on your network:

```bash
python test_flask.py --host 0.0.0.0 --port 8080
```

**Warning:** This makes the server accessible to other machines. Use with caution.

### Combined Options

You can combine multiple options:

```bash
python test_flask.py --path /home/user/documents --port 8080 --host 0.0.0.0
```

### Command-line Help

For more information on available options:

```bash
python test_flask.py --help
```

## Configuration

You can modify the default settings by editing the variables at the top of `test_flask.py`:

```python
PATH_TO_FILES = '~/ari/data/'  # Default path to explore
DEFAULT_PORT = 5000             # Default port number
```

## Project Structure

```
test_flask/
├── README.md              # This file
├── requirements.txt       # Python dependencies
├── test_flask.py         # Main Flask application
└── templates/
    └── index.html        # HTML template for the file explorer
```

## Security Notes

- The application includes a security check to prevent accessing files outside the specified root directory
- Files and directories that cannot be accessed due to permissions are automatically skipped
- When hosting with `--host 0.0.0.0`, ensure you're on a trusted network

## Troubleshooting

### Port Already in Use

If you see an error about the port being in use, try a different port:

```bash
python test_flask.py --port 8080
```

### Permission Denied Errors

If you encounter permission errors when browsing directories:
- The application will skip files/directories it cannot access
- Make sure you have read permissions for the root directory you're trying to explore

### Module Not Found: flask

If you get `ModuleNotFoundError: No module named 'flask'`, make sure you've installed the requirements:

```bash
pip install -r requirements.txt
```

## Stopping the Server

To stop the Flask server, press `Ctrl+C` in the terminal where it's running.

## Development

To run in debug mode (automatically reloads on code changes), the server is already configured with `debug=True`. Simply run the application as normal.

## Requirements

- Python 3.6 or higher
- Flask 2.0.0 or higher

## License

See the main APERO-utils repository for license information.

## Author

Neil Cook (neil.cook@umontreal.ca)

Created: January 20, 2026
