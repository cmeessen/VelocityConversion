from setuptools import setup
from setuptools import find_packages
from pkg_resources import resource_filename
from VelocityConversion import __version__ as VERSION
import versioneer

VERSION = versioneer.get_version()

# METADATA
NAME = 'velocityconversion'
MODULE = 'VelocityConversion'
AUTHOR = 'Christian Meeßen'
AUTHOR_EMAIL = 'christian.meessen@gfz-potsdam.de'
URL = 'https://github.com/cmeessen/VelocityConversion'
DESCRIPTION = 'Conversion of seismic velocities to temperature and density'
try:
    with open(resource_filename(MODULE, '../README.md'), 'r') as fh:
        LONG_DESCRIPTION = fh.read()
except ImportError:
    with open('README.md') as fh:
        LONG_DESCRIPTION = fh.read()
LONG_DESCRIPTION_TYPE = 'text/markdown'

PACKAGES = [MODULE]
PACKAGE_DIR = {MODULE: MODULE}
PACKAGE_DATA = {MODULE: ['*.csv']}

CLASSIFIERS = [
    'Natural Language :: English',
    'Programming Language :: Python :: 3',
    'Programming Language :: Python :: 3.7',
    'Programming Language :: Python :: 3.8',
    'Programming Language :: Python :: 3.9',
    'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
    'Operating System :: OS Independent',
    'Topic :: Scientific/Engineering :: Physics',
]

# DEPENDENCIES
INSTALL_REQUIRES = [
    'numpy',
]

if __name__ == '__main__':
    setup(
        name=NAME,
        version=VERSION,
        author=AUTHOR,
        author_email=AUTHOR_EMAIL,
        description=DESCRIPTION,
        long_description=LONG_DESCRIPTION,
        long_description_content_type=LONG_DESCRIPTION_TYPE,
        url=URL,
        packages=PACKAGES,
        package_dir=PACKAGE_DIR,
        package_data=PACKAGE_DATA,
        use_package_data=True,
        classifiers=CLASSIFIERS,
        install_requires=INSTALL_REQUIRES,
    )
