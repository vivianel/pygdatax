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
    'pyFAI~=0.20.0',
    'nexusformat~=0.6.1',
    'Ipython',
    'openpyxl~=3.0.7'
]

required_packages = [
    'silx',
    'matplotlib',
    'scipy',
    'pyFAI',
    'nexusformat',
    'openpyxl'
]

setup(
    name='pygdatax',
    version='0.0.9',
    url='',
    license='',
    author='achennev',
    author_email='alexis.chenneviere@cea.fr',
    description='Python Based Progamme for SAXS data treatment on XEUSS diffractometer',
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    python_requires='>=3.6',
    install_requires=required_packages,
    include_package_data=True,
    entry_points={'console_scripts': ['pygdatax_gui=pygdatax.gui_survey:main']}
)
# packages=[find_packages(where="src"])
# package_dir = {"": "src"},
