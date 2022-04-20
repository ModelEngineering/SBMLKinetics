'''Python Abstraction for a Reaction.'''

from SBMLKinetics.common import constants as cn
from SBMLKinetics.common.kinetic_law import KineticLaw
from SBMLKinetics.common import util

import numpy as np
import libsbml


class Reaction(object):

  def __init__(self, libsbml_reaction, function_definitions=None):
    """
    :param libsbml.Reaction libsbml_reaction:
    :param function_definitions list-FunctionDefinition:
    """
    self.reaction = libsbml_reaction
    # List of species reference
    self.reactants =  [self.reaction.getReactant(n)
        for n in range(self.reaction.getNumReactants())]
    # List of species reference
    self.products =  [self.reaction.getProduct(n)
        for n in range(self.reaction.getNumProducts())]
    self.kinetic_law = KineticLaw(
        #self.reaction.getKineticLaw(), self, function_definitions=None)
        self.reaction.getKineticLaw(), self, function_definitions=function_definitions)
    self.id = self.reaction.getId()

  def getId(self):
    return self.id

  def __repr__(self):
    reactant_str = " + ".join([r.getSpecies() for r in self.reactants])
    product_str = " + ".join([p.getSpecies() for p in self.products])
    if self.kinetic_law.expanded_formula is not None:
      kinetic_str = self.kinetic_law.expanded_formula
    else:
      kinetic_str = self.kinetic_law.formula
    return "%s -> %s; %s" % (reactant_str, product_str, kinetic_str)
