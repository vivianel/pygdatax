from setuptools import setup, find_packages

with open('requirements.txt', 'r') as fid:
    text = fid.readline()
    required_packages = text.splitlines()


setup(
    name='pydatax',
    version='0.0.1',
    url='',
    license='',
    author='achennev',
    author_email='alexis.chenneviere@cea.fr',
    description='Python Based Progamme for SAXS data treatment on XEUSS diffractometer',
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    install_requires=required_packages
)
