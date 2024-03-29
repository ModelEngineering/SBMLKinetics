"""
Tests for Reactions
"""
from SBMLKinetics.common import constants as cn
from SBMLKinetics.common.simple_sbml import SimpleSBML
from SBMLKinetics.common import simple_sbml
from SBMLKinetics.common.function_definition import FunctionDefinition
from tests.common import helpers

import copy
import libsbml
import numpy as np
import unittest


IGNORE_TEST = False
IS_PLOT = False


#############################
# Tests
#############################
class TestFunctionDefinition(unittest.TestCase):

  def setUp(self):
    self.simple = helpers.getSimple_BIOMD56()
    self.function_definition = FunctionDefinition(
        self.simple.model.getFunctionDefinition(0))

  def testConstructor(self):
    if IGNORE_TEST:
      return
    self.assertEqual(len(self.function_definition.argument_names), 4)


if __name__ == '__main__':
  unittest.main()
