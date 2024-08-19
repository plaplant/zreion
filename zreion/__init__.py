"""Top-level module for zreion"""

from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("zreion")
except PackageNotFoundError:
    __version__ = "unknown"

__author__ = """Paul La Plante"""
__email__ = "paul.laplante@unlv.edu"
