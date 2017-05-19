from setuptools import setup

setup(
    name='cabsDock',
    version='0.99',
    packages=['cabsDock'],
    url='https://github.com/mkurc/pyCabsDock',
    license='GNU GENERAL PUBLIC LICENSE v.3',
    author='Mateusz Kurcinski',
    author_email='mkurc@chem.uw.edu.pl',
    description='CABSDock in python',
    requires=['numpy'],
    entry_points={'console_scripts': ['cabsDock = cabsDock.__main__:run_job']},
    package_data={'cabsDock': ['data/*.dat']}
)
