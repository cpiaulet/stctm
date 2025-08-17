from setuptools import setup, find_packages
import os

def read_version():
    with open(os.path.join("stctm", "__init__.py")) as f:
        for line in f:
            if line.startswith("__version__"):
                delim = '"' if '"' in line else "'"
                return line.split(delim)[1]

setup(
    name="stctm",
    version=read_version(),
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        # put your dependencies here, e.g.:
        # "numpy", "matplotlib", "astropy", "emcee",
    ],
    entry_points={
        "console_scripts": [
            "stctm_TLS=stctm.cli.stctm_tls:main",
            # later you’ll add exotune too:
            # "stctm_exotune=stctm.cli.stctm_exotune:main",
        ],
    },
)
