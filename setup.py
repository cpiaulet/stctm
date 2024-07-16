from setuptools import setup
from pip.req import parse_requirements

with open("stctm/version.py", "r") as f:
    exec(f.read())

setup(name='stctm',
      version=__version__,
      description='Stellar contamination retrievals and modeling on small planet transmission spectra',
      url='http://github.com/cpiaulet/stctm',
      author='Caroline Piaulet-Ghorayeb',
      author_email='caroline.piaulet@umontreal.ca',
      license='GNU GPL v3.0',
      packages=['stctm'],
      install_requires=parse_requirements("requirements.txt")
      zip_safe=False)