# Kinetics-Validator
Repository of tests to validate kinetics expressions in SBML models.

## Setup environment
- Install conda
- Create a new environment
- ``conda env create -f environment.yml``
  - Sometimes you may need to separately use ``pip install python-libsbml``.
  - Also do ``conda install -c anaconda nose ``

## Running Codes
To use this environment, ``conda activate base``, where
base is the name. To exit the environment, ``conda deactivate``.
Now, do the following:

- Change directory to the top level folder of the repositor, ``kinetics_validator``.
- ``source setup_run.sh``

After completing these steps, you can run test codes. For example, to run tests for ``simple_sbml``, use
``python tests/common/test_simple_sbml.py``.
