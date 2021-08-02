from setuptools import setup, find_packages
# import os

# with open('./requirements.txt', 'r') as fid:
#     text = fid.readline()
#     required_packages = text.splitlines()

required_packages = [
    'numpy~=1.21.1',
    'silx',
    'fabio~=0.12.0',
    'matplotlib~=3.4.2',
    'scipy~=1.7.1',
    'pyFAI~=0.20.0',
    'nexusformat~=0.6.1',
    'PyQt5',
    'Ipython',
    'openpyxl~=3.0.7'
]

setup(
    name='pydatax',
    version='0.0.1',
    url='',
    license='',
    author='achennev',
    author_email='alexis.chenneviere@cea.fr',
    description='Python Based Progamme for SAXS data treatment on XEUSS diffractometer',
    package_dir={"": "src"},
    packages=['pygdatax'],
    python_requires='>=3.6',
    install_requires = required_packages,
    include_package_data=True
)
# packages=[find_packages(where="src"])
# package_dir = {"": "src"},
