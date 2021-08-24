"""
Tests for kinetics_classification.py
"""
from src import kinetics_classification
from src.common import constants as cn
from src.common import kinetic_law
from src.common.simple_sbml import SimpleSBML
from src.common.kinetic_law import KineticLaw
from tests.common import helpers
from sympy import *

import copy
import libsbml
import numpy as np
import unittest

import sys, os 


IGNORE_TEST = False
sys.stdout = open(os.devnull, 'w') #try to block the print from the main() function

#############################
# Tests
#############################
class TestKineticsClassification(unittest.TestCase):
  def testPrint1(self):
    # Test print types_name is a list of string
    if IGNORE_TEST:
      return    
    initial_model_indx = 5
    final_model_indx = 6
    #check for biomodel6
    (types_name, _, _) = kinetics_classification.main(initial_model_indx, final_model_indx)
    test = all(isinstance(item, str) for item in types_name)
    self.assertTrue(test) 

  def testPrint2(self):
    # Test print rxn_classification_num is a list of integers
    if IGNORE_TEST:
      return    
    initial_model_indx = 5
    final_model_indx = 6
    #check for biomodel6
    (_, rxn_classification_num, _) = kinetics_classification.main(initial_model_indx, final_model_indx)
    test = all(isinstance(item, int) for item in rxn_classification_num)
    self.assertTrue(test) 

  def testPrint3(self):
    # Test print rxn_num is a list of integers
    if IGNORE_TEST:
      return    
    initial_model_indx = 5
    final_model_indx = 6
    #check for biomodel6
    (_, _, rxn_num) = kinetics_classification.main(initial_model_indx, final_model_indx)
    test = isinstance(rxn_num, int)
    self.assertTrue(test) 

if __name__ == '__main__':
  unittest.main()