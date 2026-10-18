"""Package version, read from the installed distribution metadata."""

from importlib import metadata

version = metadata.version(__package__)
