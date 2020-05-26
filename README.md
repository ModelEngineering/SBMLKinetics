# Kinetics Validator
Repository of tests to validate kinetics expressions in SBML models.

## Setup environment
- Install [spyder3](http://www.psych.mcgill.ca/labs/mogillab/anaconda2/lib/python2.7/site-packages/spyder/doc/installation.html)
- Create a virtual environment for the project
  - ``cd kinetics_validator``
  - ``python -m venv kv``
  - ``source kv/bin/activate`` (use "\" in windows)
  - ``pip install -r requirements.txt``
  - ``deactivate``

To verify the setup:
- ``source kv/bin/activate`` (use "\" in windows)
- ``export PYTHONPATH=`pwd``
- ``python tests/common/test_simple_sbml.py``. The
tests should run without error.

## Running Codes
- ``cd kinetics_validator``
- ``source kv/bin/activate`` (use "\" in windows)
When you're done, use ``deactivate``.
