from setuptools import setup

__version__ = '0.9.3'

setup(
    name='CABS',
    version=__version__,
    packages=['CABS'],
    url='https://bitbucket.org/lcbio/cabsdock',
    license='free for non-commercial users',
    author='Laboratory of Computational Biology',
    author_email='mkurc@cnbc.uw.edu.pl',
    install_requires=['numpy', 'matplotlib>=2.0', 'requests'],
    description='CABS in python',
    entry_points={
        'console_scripts': [
            'cabsDock = CABS.__main__:run_dock',
            'cabsFlex = CABS.__main__:run_flex'
        ]
    },
    package_data={'CABS': ['data/*.dat']}
)
