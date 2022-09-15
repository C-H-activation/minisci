"""Bash helper functions."""

try:
    import os
except ImportError:
    raise ImportError("Module os is required. Please install module os.")
try:
    import shutil
except ImportError:
    raise ImportError("Module shutil is required. Please install module shutil.")


def copy_file(source: str, destination: str):
    """Copy file from source (path) to destination (path)."""
    shutil.copy2(source, destination)


def remove_file(path: str):
    """Remove file safely."""
    if os.path.exists(path):
        os.remove(path)
