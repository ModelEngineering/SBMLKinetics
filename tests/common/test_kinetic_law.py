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


#############################
# Tests
#############################
class TestKineticLaw(unittest.TestCase):

  def setUp(self):
    self.simple = helpers.getSimple()
    self.law = self.simple.reactions[1].kinetic_law

  def testConstructor(self):
    if IGNORE_TEST:
      return
    self.assertTrue(
        isinstance(self.law, KineticLaw))


if __name__ == '__main__':
  unittest.main()
