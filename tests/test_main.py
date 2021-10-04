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
import pandas as pd
import math

IGNORE_TEST = False
#sys.stdout = open(os.devnull, 'w') #try to block the print from the main() function

#############################
# Tests
#############################
class TestKineticsClassification(unittest.TestCase):
  def testClassification1(self):
    # Test all the column names
    if IGNORE_TEST:
      return    
    initial_model_indx = 5
    final_model_indx = 6
    #check for biomodel6
    (df_classification, _, _) = kinetics_classification.main(initial_model_indx, final_model_indx)
    test = all(item in df_classification.columns for item in kinetics_classification.COLUMN_NAME_df_classification)
    self.assertTrue(test)

  def testClassification2(self):
    # Test whether there is at least one row
    if IGNORE_TEST:
      return    
    initial_model_indx = 5
    final_model_indx = 6
    #check for biomodel6
    (df_classification, _, _) = kinetics_classification.main(initial_model_indx, final_model_indx)
    self.assertTrue(len(df_classification.index)>0) 

  def testClassification3(self):
    # Test all the elements of df_classification are lists of strings
    if IGNORE_TEST:
      return    
    initial_model_indx = 5
    final_model_indx = 6
    #check for biomodel6
    (df_classification, _, _) = kinetics_classification.main(initial_model_indx, final_model_indx)
    list_classification = []
    for i in range(len(df_classification.columns)):
      list_classification += df_classification.iloc[:,i].tolist()
    test = all(isinstance(item, str) for item in list_classification)
    self.assertTrue(test)

  def testGenStat1(self):
    # Test all the column names
    if IGNORE_TEST:
      return    
    initial_model_indx = 5
    final_model_indx = 6
    #check for biomodel6
    (_, df_gen_stat, _) = kinetics_classification.main(initial_model_indx, final_model_indx)
    test = all(item in df_gen_stat.columns for item in kinetics_classification.COLUMN_NAME_df_gen_stat)
    self.assertTrue(test)

  def testGenStat2(self):
    # Test whether there is at least one row
    if IGNORE_TEST:
      return    
    initial_model_indx = 5
    final_model_indx = 6
    #check for biomodel6
    (_, df_gen_stat, _) = kinetics_classification.main(initial_model_indx, final_model_indx)
    self.assertTrue(len(df_gen_stat.index)>0) 

  def testGenStat3(self):
    # Test column 'Classification Names' of df_gen_stat is a list of string
    if IGNORE_TEST:
      return    
    initial_model_indx = 5
    final_model_indx = 6
    #check for biomodel6
    (_, df_gen_stat, _) = kinetics_classification.main(initial_model_indx, final_model_indx)
    list_gen_stat_classification = df_gen_stat[kinetics_classification.COLUMN_NAME_df_gen_stat[0]].tolist()
    test = all(isinstance(item, str) for item in list_gen_stat_classification)
    self.assertTrue(test) 

  def testGenStat4(self):
    # Test column 'Percentage' of df_gen_stat is a list of floating numbers
    if IGNORE_TEST:
      return    
    initial_model_indx = 5
    final_model_indx = 6
    #check for biomodel6
    (_, df_gen_stat, _) = kinetics_classification.main(initial_model_indx, final_model_indx)
    list_gen_stat_percentage = df_gen_stat[kinetics_classification.COLUMN_NAME_df_gen_stat[1]].tolist()
    test = all(isinstance(item, float) for item in list_gen_stat_percentage)
    self.assertTrue(test) 

  def testGenStat5(self):
    # Test column 'Percentage' of df_gen_stat does not have nan values
    if IGNORE_TEST:
      return    
    initial_model_indx = 5
    final_model_indx = 6
    #check for biomodel6
    (_, df_gen_stat, _) = kinetics_classification.main(initial_model_indx, final_model_indx)
    list_gen_stat_percentage = df_gen_stat[kinetics_classification.COLUMN_NAME_df_gen_stat[1]].tolist()
    test = any(math.isnan(item) for item in list_gen_stat_percentage)
    self.assertFalse(test)

  def testMolStat1(self):
    # Test all the column names
    if IGNORE_TEST:
      return    
    initial_model_indx = 5
    final_model_indx = 6
    #check for biomodel6
    (_, _ , df_mol_stat) = kinetics_classification.main(initial_model_indx, final_model_indx)
    test = all(item in df_mol_stat.columns for item in kinetics_classification.COLUMN_NAME_df_mol_stat)
    self.assertTrue(test)
 
  def testMolStat2(self):
    # Test whether there is at least one row
    if IGNORE_TEST:
      return    
    initial_model_indx = 5
    final_model_indx = 6
    #check for biomodel6
    (_, _, df_mol_stat) = kinetics_classification.main(initial_model_indx, final_model_indx)
    self.assertTrue(len(df_mol_stat.index)>0) 

  def testMolStat3(self):
    # Test column 'SBMLid' of df_mol_stat is a list of string
    if IGNORE_TEST:
      return    
    initial_model_indx = 5
    final_model_indx = 6
    #check for biomodel6
    (_, _, df_mol_stat) = kinetics_classification.main(initial_model_indx, final_model_indx)
    list_gen_stat_classification = df_mol_stat[kinetics_classification.COLUMN_NAME_df_mol_stat[0]].tolist()
    test = all(isinstance(item, str) for item in list_gen_stat_classification)
    self.assertTrue(test)

  def testMolStat4(self):
    # Test column 'Reaction#' of df_mol_stat is a list of integers
    if IGNORE_TEST:
      return    
    initial_model_indx = 5
    final_model_indx = 6
    #check for biomodel6
    (_, _, df_mol_stat) = kinetics_classification.main(initial_model_indx, final_model_indx)
    list_gen_stat_classification = df_mol_stat[kinetics_classification.COLUMN_NAME_df_mol_stat[1]].tolist()
    test = all(isinstance(item, int) for item in list_gen_stat_classification)
    self.assertTrue(test)

  def testMolStat5(self):
    # Test columns other than 'SBMLid' and 'Reaction#' of df_mol_stat are lists of floating numbers
    if IGNORE_TEST:
      return    
    initial_model_indx = 5
    final_model_indx = 6
    #check for biomodel6
    (_, _, df_mol_stat) = kinetics_classification.main(initial_model_indx, final_model_indx)
    list_gen_stat_others = []
    for i in range(2,len(df_mol_stat.columns)):
      list_gen_stat_others += df_mol_stat.iloc[:,i].tolist()
    test = all(isinstance(item, float) for item in list_gen_stat_others)
    self.assertTrue(test)

  def testMolStat6(self):
    # Test whether all the numbers are not nan values
    if IGNORE_TEST:
      return    
    initial_model_indx = 5
    final_model_indx = 6
    #check for biomodel6
    (_, _, df_mol_stat) = kinetics_classification.main(initial_model_indx, final_model_indx)
    list_gen_stat_others = []
    for i in range(1,len(df_mol_stat.columns)):
      list_gen_stat_others += df_mol_stat.iloc[:,i].tolist()
    test = any(math.isnan(item) for item in list_gen_stat_others)
    self.assertFalse(test)

  def testMolStat7(self):
    # Test the sum of percentage of all types for each model in df_mol_stat is always no less than 1
    if IGNORE_TEST:
      return    
    initial_model_indx = 5
    final_model_indx = 6
    #check for biomodel6
    (_, _, df_mol_stat) = kinetics_classification.main(initial_model_indx, final_model_indx)
    list_gen_stat_others = []
    flag = 1
    for i in range(len(df_mol_stat)):
      sum = 0
      list_gen_stat_others += df_mol_stat.iloc[i,2:len(df_mol_stat.columns)].tolist()
      for j in range(len(list_gen_stat_others)):
        sum += list_gen_stat_others[j]
      if sum < 1:
        flag = 0
    self.assertTrue(flag)


if __name__ == '__main__':
  unittest.main()
