from setuptools import setup, find_packages
import os.path
import sys
import codecs

with open("README.md", "r") as fh:
    long_description = fh.read()

# The following two methods were copied from
# https://packaging.python.org/guides/single-sourcing-package-version/#single-sourcing-the-version
def read(rel_path):
    here = os.path.abspath(os.path.dirname(__file__))
    with codecs.open(os.path.join(here, rel_path), 'r') as fp:
        return fp.read()

def get_version(rel_path):
    for line in read(rel_path).splitlines():
        if line.startswith('__version__'):
            print (line)
            delim = '"' if '"' in line else "'"
            print ('delim = ', delim)
            return line.split(delim)[1]
    else:
        raise RuntimeError("Unable to find version string.")

setup(
    name='SBMLKinetics',
    packages=find_packages(),
    version=get_version('SBMLKinetics/_version.py'),
    description='Analyze SBML kinetics.',
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='Jin Xu, Joseph Hellerstein',
    author_email='jxu2019@uw.edu',
    url='https://github.com/ModelEngineering/SBMLKinetics',
    license='MIT License',
    install_requires=[
        'coverage',
        'matplotlib',
        'nose',
        'numpy',
        'pandas',
        'pylint',
        'python-libsbml',
        'pip>20',
        'sympy',
        'tellurium!=2.2.7',
        'urllib3',
        'seaborn',
        'xlsxwriter',
        'openpyxl',
    ],

    #scripts=[''],# The name of your scipt, and also the command you'll be using for calling it
    scripts=[
    'SBMLKinetics/kinetics_classification.py',
    'SBMLKinetics/common/constants.py',
    'SBMLKinetics/common/exceptions.py',
    'SBMLKinetics/common/function_definition.py',
    'SBMLKinetics/common/helpers.py',
    'SBMLKinetics/common/kinetic_law.py',
    'SBMLKinetics/common/msgs.py',
    'SBMLKinetics/common/reaction.py',
    'SBMLKinetics/common/simple_sbml.py',
    'SBMLKinetics/common/util.py',
    ],
    include_package_data=True,
    data_files=[('SBMLKinetics/data', ['SBMLKinetics/data/biomodels.zip']), 
                ('SBMLKinetics/data', ['SBMLKinetics/data/curated.zip']), 
                ('SBMLKinetics/data', ['SBMLKinetics/data/signalling.zip']),
                ('SBMLKinetics/data', ['SBMLKinetics/data/metabolic.zip']), 
                ('SBMLKinetics/data', ['SBMLKinetics/data/homo_sapiens.zip']), 
                ('SBMLKinetics/data', ['SBMLKinetics/data/non_homo.zip']),
                ('SBMLKinetics/data', ['SBMLKinetics/data/cellular_organisms.zip']), 
                ('SBMLKinetics/data', ['SBMLKinetics/data/Mammalia.zip']), 
                ('SBMLkinetics/data', ['SBMLKinetics/data/Mus_musculus.zip']),
                ('SBMLkinetics/data', ['SBMLKinetics/data/Saccharomyces_cerevisiae.zip'])],
    classifiers=[
       'License :: OSI Approved :: MIT License',
       'Programming Language :: Python :: 3.8',
       'Programming Language :: Python :: 3.9',
       'Programming Language :: Python :: 3.10',
       'Operating System :: OS Independent',
    ],
    python_requires='>=3.8',
)
