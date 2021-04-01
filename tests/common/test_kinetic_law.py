"""
Tests for Kinetic Law
"""
from src.common import constants as cn
from src.common.simple_sbml import SimpleSBML
from src.common.kinetic_law import KineticLaw
from tests.common import helpers

import copy
import libsbml
import numpy as np
import unittest


IGNORE_TEST = False
IS_PLOT = False
NUM_LAW = 3


class MockFunctionDefinition():
  """Used to mock function defintions."""

  def __init__(self, id, argument_names, body):
    self.id = id
    self.argument_names = argument_names
    self.body = body

  def __repr__(self):
    argument_call = ",".join(self.argument_names)
    call_str = "%s(%s)" % (self.id, argument_call)
    return "%s: %s" % (call_str, self.body)


#############################
# Tests
#############################
class TestKineticLaw(unittest.TestCase):

  def setUp(self):
    self.simple = helpers.getSimple_BIOMD6()
    self.laws = [self.simple.reactions[i].kinetic_law for i in range(NUM_LAW)]

  def testConstructor(self):
    if IGNORE_TEST:
      return
    self.assertTrue(
         isinstance(self.laws[1], KineticLaw))
    trues = [isinstance(s, str) for s in self.laws[1].symbols]
    self.assertTrue(all(trues))

  def testSymbol(self):
    if IGNORE_TEST:
      return
    def checkSubset(subset,superset):
      true = set(subset).issubset(set(superset))
      self.assertTrue(true)
    # const only
    checkSubset(['kappa'],self.laws[0].symbols)
    # simple kinetic law k*A
    checkSubset(['k6','u'],self.laws[1].symbols)
    # complex kinetic law with function pow()
    checkSubset(['z','u','k4','k4prime'],self.laws[2].symbols)

  def testExpandFormula1(self):
    # Replacement for SBML reactions
    if IGNORE_TEST:
      return
    simple = helpers.getSimple_BIOMD56()
    kinetic_law = None
    for fd in simple.function_definitions:
      for reaction in simple.reactions:
        if fd.id in reaction.kinetic_law.formula:
          kinetic_law = KineticLaw(reaction.kinetic_law.libsbml_kinetics,
              reaction, function_definitions=simple.function_definitions)
          break
      if kinetic_law is not None:
        break
    if kinetic_law is None:
      raise RuntimeError("Did not find an embedded function.")
    self.assertIsNotNone(kinetic_law.expanded_formula)
    kinetic_law_arguments = ["Vaiep", "Jaiep", "1", "IE"]
    kinetic_law.expandFormula(simple.function_definitions)
    for argument in kinetic_law_arguments:
      self.assertTrue(argument in kinetic_law.expanded_formula)

  def mkKineticLawWithFormula(self, formula):
    kinetic_law = KineticLaw(self.laws[0].libsbml_kinetics, None)
    kinetic_law.formula = formula
    return kinetic_law

  def testExpandFormula2(self):
    # Simple replacement
    if IGNORE_TEST:
      return
    function_definitions = [
        MockFunctionDefinition("aa", ["x", "y"], "x + y")
        ]
    kinetic_law = self.mkKineticLawWithFormula("aa(1, 2)")
    kinetic_law.expandFormula(function_definitions)
    self.assertEqual(kinetic_law.expanded_formula, "1 + 2")

  def testExpandFormula4(self):
    # Test replacement with other embedded expressions
    if IGNORE_TEST:
      return
    function_definitions = [
        MockFunctionDefinition("aa", ["x", "y"], "x + y")
        ]
    kinetic_law = self.mkKineticLawWithFormula("exp(4) + aa(1, 2) + sin(3)")
    kinetic_law.expandFormula(function_definitions)
    self.assertEqual(kinetic_law.expanded_formula, "exp(4) + 1 + 2 + sin(3)")

  def testExpandFormula4(self):
    # Test replacing multiple functions
    if IGNORE_TEST:
      return
    function_definitions = [
        MockFunctionDefinition("aa", ["x", "y"], "x + y"),
        MockFunctionDefinition("bb", ["x", "y"], "x*y")
        ]
    kinetic_law = self.mkKineticLawWithFormula("2 + aa(1, 2) + bb(x, z)")
    kinetic_law.expandFormula(function_definitions)
    self.assertEqual(kinetic_law.expanded_formula, "2 + 1 + 2 + x*z")

  def testExpandFormula5(self):
    # Test recursive replacements
    if IGNORE_TEST:
      return
    function_definitions = [
        MockFunctionDefinition("aa", ["x", "y"], "x + y"),
        MockFunctionDefinition("bb", ["x", "y"], "bb + aa(x, y)"),
        MockFunctionDefinition("cc", ["x", "y"], "cc + bb(x, y)")
        ]
    kinetic_law = self.mkKineticLawWithFormula("kl + cc(1, 2)")
    kinetic_law.expandFormula(function_definitions)
    self.assertEqual(kinetic_law.expanded_formula, "kl + cc + bb + 1 + 2")


if __name__ == '__main__':
  unittest.main()

