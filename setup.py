from setuptools import setup

setup(
    name='CABS',
    version='1.0',
    packages=['cabsDock'],
    url='https://bitbucket.org/lcbio/cabsdock',
    license='free for non-commercial users',
    author='Laboratory of Computational Biology',
    author_email='mkurc@cnbc.uw.edu.pl',
    description='CABS in python',
    install_requires=['numpy', 'matplotlib>=2.0'],
    entry_points={'console_scripts': [  'cabsDock = cabsDock.__main__:run_dock',
                                        'cabsFlex = cabsDock.__main__:run_flex']},
    package_data={'cabsDock': ['data/*.dat']}
)
