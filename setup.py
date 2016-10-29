"""
Good place to read about packaging:
https://python-packaging-user-guide.readthedocs.org/en/latest/distributing.html
"""

from setuptools import setup, find_packages

setup(name='ODIN',
      version='1.0',
      packages=find_packages(),
      install_requires=['numpy==1.11.1']
      )