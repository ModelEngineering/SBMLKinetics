"""
This is a tutorial on using libsbml and the python
wrappers.
Make sure that you have setup your PYTHONPATH environment
variable as described in the github repository.
"""

# Import the required files
from src.common.simple_sbml import SimpleSBML
import src.common.simple_sbml as simple_sbml
import src.common.constants as cn

import numpy as np
import os

# Create an SBML model. We'll use the model
# tests/common/test_file.xml
path = os.path.join(cn.TEST_DIR, "common")
path = os.path.join(path, "test_file.xml")
simple = SimpleSBML(path)  # Creates a model

# The model has instance variables for the following:
#  simple.species - list of libsbml Species objects
#  simple.reactions - list of Reaction objects as
#                     described in src/common/reaction.py
#  simple.reactions - list of Reaction objects as
#                     described in src/common/reaction.py
#  simple.parameters- list of libsbml Parameters objects

# Here's an example of using these objects
names = ""
for spc in simple.species:
  names = names + " " + spc.getId()
print("\n***Species names***\n  %s" % names)

# Reactions have a more complex structure.
#  reactants - list of libsbl SpeciesReference
#  products - list of libsbl SpeciesReference
#  kinetic_law - KinetcLaw
# KineticLaw object that provides
#  formula - string of kinetics formula
#  symbols - list of chemical species and parameters
print("\n***Reactions w/o stoichiometries.")
for reaction in simple.reactions:
  reactant_stg = " + ".join(
      [r.getSpecies() for r in reaction.reactants])
  product_stg = " + ".join(
      [p.getSpecies() for p in reaction.products])
  print("%s -> %s; %s" % (
      reactant_stg, product_stg,
      reaction.kinetic_law.formula))

# The code below shows how to iterate through all
# of BioModels. The iterator returns an object
# that contains the file name and a SimpleSBML object.
print("\n***Iterating through BioModels.")
iterator = simple_sbml.modelIterator(initial=0, final=10)
for item in iterator:
    print("File %s has %d reactions." % (item.filename, len(item.model.reactions)))
