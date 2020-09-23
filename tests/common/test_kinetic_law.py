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
    #self.simple = helpers.getSimple()
    self.simple = helpers.getSimple_BIOMD6()
    global law 
    law = []
    for i in range(3):
      law.append(self.simple.reactions[i].kinetic_law)

  def testConstructor(self):
    if IGNORE_TEST:
      return
    self.assertTrue(
         isinstance(law[1], KineticLaw))
    trues = [isinstance(s, str) for s in law[1].symbols]
    self.assertTrue(all(trues))

  def testSymbol(self):
    def checkSubset(subset,superset):
      true = set(subset).issubset(set(superset))
      self.assertTrue(true)

    # const only
    checkSubset(['kappa'],law[0].symbols)

    # simple kinetic law k*A
    checkSubset(['k6','u'],law[1].symbols)

    # complex kinetic law with function pow()
    checkSubset(['z','u','k4','k4prime'],law[2].symbols)


if __name__ == '__main__':
  unittest.main()

