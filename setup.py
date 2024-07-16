from setuptools import setup

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
      install_requires=['numpy', 'scipy', 'emcee', 'corner', 'astropy'],
      zip_safe=False)