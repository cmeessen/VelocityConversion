from setuptools import setup
from setuptools import find_packages
from pkg_resources import resource_filename


# METADATA
NAME = 'velocityconversion-cmeessen'
MODULE = 'VelocityConversion'
VERSION = '1.1.0'
AUTHOR = 'Christian Meeßen'
AUTHOR_EMAIL = 'christian.meessen@gfz-potsdam.de'
MAINTAINER = 'Christian Meeßen'
MAINTAINER_EMAIL = 'christian.meessen@gfz-potsdam.de'
URL = 'https://github.com/cmeessen/VelocityConversion'
DESCRIPTION = 'Conversion of seismic velocities to temperature and density'
try:
    with open(resource_filename(MODULE, '../README.md'), 'r') as fh:
        LONG_DESCRIPTION = fh.read()
except ImportError:
    with open('README.md') as fh:
        LONG_DESCRIPTION = fh.read()
LONG_DESCRIPTION_TYPE = 'text/markdown'
PACKAGE_DATA = find_packages()
CLASSIFIERS = [
    'Natural Language :: English',
    'Programming Language :: Python :: 2'
    'Programming Language :: Python :: 2.7',
    'Programming Language :: Python :: 3',
    'Programming Language :: Python :: 3.7',
    'License :: OSI Approved :: GNU GPL-3.0',
    'Operating System :: OS Independent',
    'Topic :: Geophysics',
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
        maintainer=MAINTAINER,
        maintainer_email=MAINTAINER_EMAIL,
        description=DESCRIPTION,
        long_description=LONG_DESCRIPTION,
        long_description_content_type=LONG_DESCRIPTION_TYPE,
        url=URL,
        packages=PACKAGE_DATA,
        classifiers=CLASSIFIERS,
        install_requires=INSTALL_REQUIRES,
    )
