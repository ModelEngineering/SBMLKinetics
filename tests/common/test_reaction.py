"""
Tests for Reactions
"""
from src.common import constants as cn
from src.common.simple_sbml import SimpleSBML
from src.common.reaction import Reaction
from tests.common import helpers

import copy
import libsbml
import numpy as np
import unittest


IGNORE_TEST = True
IS_PLOT = True


#############################
# Tests
#############################
class TestReaction(unittest.TestCase):

  def setUp(self):
    self.simple = helpers.getSimple()
    self.reactions = self.simple.reactions
    self.reaction = self.reactions[2]

  def testConstructor(self):
    if IGNORE_TEST:
      return
    def test(a_list, a_type):
      self.assertGreater(len(a_list), 0)
      self.assertTrue(isinstance(a_list[0], a_type))
    #
    test(self.reaction.reactants,
        libsbml.SpeciesReference)
    test(self.reaction.products,
        libsbml.SpeciesReference)

  def testRepr(self):
    # TESTING
    reaction_str = str(self.reaction)
    for item in ["T2R", "->", "R", "k1c", "*"]:
      self.assertTrue(item in reaction_str)


if __name__ == '__main__':
  unittest.main()
