'''Provides Information on SBML Kinetics Laws'''


from src.common import constants as cn
from src.common import util
from src.common import exceptions
from src.common import msgs

import collections
import numpy as np
import re


class KineticLaw(object):

  def __init__(self, libsbml_kinetics, reaction):
    """
    :param libsbml.KineticLaw libsbml_kinetics:
    """
    # libsbml object for kinetics
    self.libsbml_kinetics = libsbml_kinetics
    # String version of chemical formula
    self.formula = self.libsbml_kinetics.getFormula()
    # Parameters and chemical species
    # FIXME: Why doesn't this work
    try:
        self.symbols = self._getSymbols()
    except Exception:
        msgs.warn("Could not acquire symbols.")
        self.symbols = None
    # Reaction for the kinetics law
    self.reaction = reaction

  def __repr__(self):
    return self.formula

  def expandFormula(self, function_definitions):
    """
    Expands the kinetics formula, replacing function definitions
    with their body.

    Parameters
    ----------
    function_definitions: list-FunctionDefinition
    
    Returns
    -------
    str
    """
    expansion = str(self.formula)
    for fd in function_definitions:
      # Find the function calls
      calls = re.findall(r'{}\(.*?\)'.format(fd.id), self.formula)
      for call in calls:
        # Find argument call. Ex: '(a, b)'
        call_arguments = re.findall(r'\(.*?\)', call)[0]
        call_arguments = call_arguments.strip()
        call_arguments = call_arguments[1:-1]  # Eliminate parentheses
        arguments = call_arguments.split(',')
        body = str(fd.body)
        for formal_arg, call_arg in zip(fd.argument_names, arguments):
          body = body.replace(formal_arg, call_arg)
        expansion = expansion.replace(call, body)
    return expansion
 

  def _getSymbols(self):
    """
    Finds the parameters and species names for the
    kinetics law. Exposing this information requires
    a recursive search of the parse tree for the
    kinetics expression.
    :return list-str:
    """
    global cur_depth
    MAX_DEPTH = 20
    cur_depth = 0
    def augment(ast_node, result):
      global cur_depth
      cur_depth += 1
      if cur_depth > MAX_DEPTH:
        raise exceptions.BadKineticsMath(self.reaction.id)
      for idx in range(ast_node.getNumChildren()):
        child_node = ast_node.getChild(idx)
        if child_node.getName() is None:
          additions = augment(child_node, result)
          result.extend(additions)
        else:
          if child_node.isFunction():
            additions = augment(child_node, result)
            result.extend(additions)
          else:
            result.append(child_node.getName())
      return result

    ast_node = self.libsbml_kinetics.getMath()
    if ast_node.getName() is None:
      result = []
    else:
      result = [ast_node.getName()]
    return augment(ast_node, result)
