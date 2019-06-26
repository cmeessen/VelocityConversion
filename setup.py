from setuptools import setup
from setuptools import find_packages

# METADATA
NAME = 'velocityconversion-cmeessen'
VERSION = '1.1.0'
AUTHOR = 'Christian Meeßen'
AUTHOR_EMAIL = 'christian.meessen@gfz-potsdam.de'
MAINTAINER = 'Christian Meeßen'
MAINTAINER_EMAIL = 'christian.meessen@gfz-potsdam.de'
URL = 'https://github.com/cmeessen/VelocityConversion'
DESCRIPTION = 'Conversion of seismic velocities to temperature and density'
with open('README.md', 'r') as fh:
    LONG_DESCRIPTION = fh.read()
LONG_DESCRIPTION_TYPE = 'text/markdown'
PACKAGE_DATA = find_packages()
CLASSIFIERS = [
    'Programming Language :: Python :: 3',
    'License :: OSI Approved :: GNU GPL-3.0',
    'Operating System :: OS Independent',
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
