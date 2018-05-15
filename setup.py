from setuptools import setup
import os.path

# Get the version:
version = {}
with open(os.path.join(os.path.abspath(os.path.dirname(__file__)), 'version.py')) as f: exec(f.read(), version)

setup(
    name = 'heterogenesis',
    version = version['__version__'],
    description = 'Description needed.',
    author = 'Georgette Tanner',
    author_email = 'medgnt@leeds.ac.uk',
    url = 'https://github.com/GeorgetteTanner/heterogenesis',
    classifiers = [
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)',
        'Programming Language :: Python :: 3'
    ],
    py_modules = ['heterogenesis_vargen','heterogenesis_varincorp','freqcalc','version'],
    install_requires = [
    'numpy>=1.12.0'
    ],
    python_requires = '>=3',
    entry_points = {
        'console_scripts': [
            'heterogenesis_vargen=heterogenesis_vargen:main',
            'heterogenesis_varincorp=heterogenesis_varincorp:main',
            'freqcalc=freqcalc:main'
        ]
    }
)
