from setuptools import setup, find_packages
import os


with open('requirements.txt') as f:
    required = f.read().splitlines()

setup(
    name='ZSeeker',  # Package name
    version='1.3',  # Version number
    description='ZSeeker is a cli tool to find the propensity of B-DNA to form Z-DNA structures. ',  # Short description
    long_description=open('README.md').read(),  # Full description from README.md
    long_description_content_type='text/markdown',  # Content type for long description
    author='Nikol Chantzi, Patsakis Michail , Provatas Kimon, Ilias Georgakopoulos Soares, Ioannis Mouratidis',  # Your name
    author_email='kap6605@psu.edu , mpp5977@psu.edu',  # Your email
    url='https://github.com/Georgakopoulos-Soares-lab/ZSeeker',  # URL to your project (optional)
    packages=find_packages(),  # Automatically find all packages in the project
    entry_points={
        'console_scripts': [
            'ZSeeker=zseeker.zseeker:main',  # Register CLI command
        ],
    },
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.10',  # Specify Python version requirement
    install_requires=[
        required
    ],
)
