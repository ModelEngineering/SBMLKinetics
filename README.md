<img src="https://api.travis-ci.org/ModelEngineering/kinetics_validator.svg?branch=master" width="100"/>

[![Coverage](https://codecov.io/gh/ModelEngineering/kinetics_validator/branch/master/graph/badge.svg)](https://codecov.io/gh/ModelEngineering/kinetics_validator)

# SBMLKinetics
SBMLKinetics is a package of tests to analyze kinetics expressions in SBML models.
If you are using any of the code, please cite the PYPI web page 
(https://pypi.org/project/SBMLKinetics/). 

## For users
### Installation

``pip install SBMLKinetics``

### Documentation
Please see the documentation at https://modelengineering.github.io/kinetics_validator/ for details.


## For developers
### Setup environment
- Install [spyder3](http://www.psych.mcgill.ca/labs/mogillab/anaconda2/lib/python2.7/site-packages/spyder/doc/installation.html)
- Clone the ``kinetics_validator`` repository using ``git clone https://github.com/ModelEngineering/kinetics_validator``
- Create a virtual environment for the project.
  - ``cd kinetics_validator``
  - ``python -m venv kv``
  - ``source kv/Scripts/activate``
(Use "\\" in windows.)
  - ``pip install -r requirements.txt``
  - ``deactivate``

To verify the setup:
- Return to the ``kinetics_evaluator`` directory.
- ``source kv/Scripts/activate``
(Use "\\" in windows.)
- ``export PYTHONPATH=`pwd` ``
- ``python tests/test_simple_sbml.py``. The
tests should run without error.
(Use "\\" in windows.)

### Running Codes
- ``cd kinetics_validator``
- ``source kv/bin/activate``
(Use "\\" in windows.)
When you're done, use ``deactivate``.

### Documentation
- ``examples/tutorial.py`` has code illustrating usage
- ``SBMLKinetics/common/*.py`` has codes for the 
SmpleSBML (``simple_sbml.py``),
Reaction (``reaction.py``),
and KineticLaw (``kinetic_law.py``).
