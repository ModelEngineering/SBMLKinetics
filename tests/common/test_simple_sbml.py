"""
Tests for simple_sbml
"""
from src.common import constants as cn
from src.common import simple_sbml
from src.common.simple_sbml import SimpleSBML
from src.common.reaction import Reaction
from src.common import util
from tests.common import helpers

import copy
import numpy as np
import os
import libsbml
import unittest
import zipfile


IGNORE_TEST = False
NO_NAME = "dummy"


#############################
# Tests
#############################
class TestSimpleSBML(unittest.TestCase):

  def setUp(self):
    self.simple = helpers.getSimple()

  def testConstructor(self):
    if IGNORE_TEST:
      return
    def test(a_list, a_type):
      self.assertGreater(len(a_list), 0)
      self.assertTrue(isinstance(a_list[0], a_type))
    #
    test(self.simple.reactions, Reaction)
    test(self.simple.species, libsbml.Species)
    test(self.simple.parameters, libsbml.Parameter)
    self.assertTrue(isinstance(self.simple.model,
        libsbml.Model))

  def testGet(self):
    if IGNORE_TEST:
      return
    def test(func, a_list):
      this_id = a_list[0].getId()
      an_object = func(this_id)
      self.assertEqual(an_object, a_list[0])
    #
    test(self.simple.getReaction, self.simple.reactions)
    test(self.simple.getSpecies, self.simple.species)
    test(self.simple.getParameter, self.simple.parameters)


class TestFunctions(unittest.TestCase):

  def testReadURL(self):
    pass

  def _testIterator(self, itr):
    for item in itr:
      self.assertTrue(isinstance(item.model,
          SimpleSBML))
    COUNT = 5
    itr = simple_sbml.modelIterator(final=COUNT)
    item_number = -1
    for item in itr:
      self.assertTrue(isinstance(item.filename, str))
      item_number = item.number
    self.assertEqual(item_number, COUNT - 1)

  def testModelIterator1(self):
    if IGNORE_TEST:
      return
    self._testIterator(simple_sbml.modelIterator(final=1))

  def testGetZipfilePath(self):
    if IGNORE_TEST:
      return
    ffiles, zipper = simple_sbml.getZipfilePaths()
    for ffile in ffiles:
      try:
        fid = zipper.open(ffile)
        fid.close()
      except:
        assertTrue(False)
    

if __name__ == '__main__':
  unittest.main()
