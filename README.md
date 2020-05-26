# Kinetics Validator
Repository of tests to validate kinetics expressions in SBML models.

## Setup environment
- Install [spyder3](http://www.psych.mcgill.ca/labs/mogillab/anaconda2/lib/python2.7/site-packages/spyder/doc/installation.html)
- Clone the ``kinetics_validator`` repository using ``git clone https://github.com/ModelEngineering/kinetics_validator``
- Create a virtual environment for the project.
  - ``cd kinetics_validator``
  - ``python -m venv kv``
  - ``source kv/bin/activate``
(Use "\\" in windows.)
  - ``pip install -r requirements.txt``
  - ``deactivate``

To verify the setup:
- Return to the ``kinetics_evaluator`` directory.
- ``source kv/bin/activate``
(Use "\\" in windows.)
- ``export PYTHONPATH=`pwd` ``
- ``python tests/common/test_simple_sbml.py``. The
tests should run without error.
(Use "\\" in windows.)

## Running Codes
- ``cd kinetics_validator``
- ``source kv/bin/activate``
(Use "\\" in windows.)
When you're done, use ``deactivate``.
