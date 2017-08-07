from setuptools import setup

setup(
    name='cabsDock',
    version='1.0',
    packages=['cabsDock'],
    url='https://bitbucket.org/lcbio/cabsdock',
    license='free for non-commercial users',
    author='Laboratory of Computational Biology',
    author_email='mkurc@cnbc.uw.edu.pl',
    description='CABSDock in python',
    requires=['numpy', 'matplotlib'],
    entry_points={'console_scripts': ['cabsDock = cabsDock.__main__:run_job']},
    package_data={'cabsDock': ['data/*.dat']}
)
