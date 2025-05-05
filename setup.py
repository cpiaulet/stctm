from setuptools import setup
import os

def read_version():
    with open(os.path.join("stctm", "__init__.py")) as f:
        for line in f:
            if line.startswith("__version__"):
                delim = '"' if '"' in line else "'"
                return line.split(delim)[1]

setup(
    # other args...
    version=read_version(),
)