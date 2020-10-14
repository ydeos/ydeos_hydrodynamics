#!/usr/bin/env python
# coding: utf-8

r"""ydeos_hydrodynamics's setup.py"""

from distutils.core import setup

import ydeos_hydrodynamics


setup(name=ydeos_hydrodynamics.__project_name__,
      version=ydeos_hydrodynamics.__version__,
      description=ydeos_hydrodynamics.__description__,
      long_description='Hydrodynamics computations: viscous, Delft, MIT, '
                       'Savitski, Holtrop-Mennen, sideforce, righting moments',
      url=ydeos_hydrodynamics.__url__,
      download_url=ydeos_hydrodynamics.__download_url__,
      author=ydeos_hydrodynamics.__author__,
      author_email=ydeos_hydrodynamics.__author_email__,
      license=ydeos_hydrodynamics.__license__,
      classifiers=[
          'Development Status :: 3 - Alpha',
          'Intended Audience :: Developers',
          'Topic :: Software Development',
          'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
          'Programming Language :: Python :: 3.7'],
      keywords='hydrodynamics',
      packages=['ydeos_hydrodynamics'])
