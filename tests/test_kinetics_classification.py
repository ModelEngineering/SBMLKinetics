"""
Tests for kinetics_classification.py
"""
from SBMLKinetics import kinetics_classification
from SBMLKinetics.common import constants as cn
from SBMLKinetics.common import kinetic_law
from SBMLKinetics.common.simple_sbml import SimpleSBML
from SBMLKinetics.common.kinetic_law import KineticLaw
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
  
  def setUp(self):
    #check for biomodel6
    initial_model_indx = 5
    final_model_indx = 6
    model_indices = range(initial_model_indx, final_model_indx+1)

    self.df_classification, self.df_gen_stat, self.df_mol_stat, self.df_gen_stat_PR, \
    self.biomodel_non_count, self.df_table_PR, self.df_table_PR_per_model \
    = kinetics_classification._dataSetToTuple(initial_model_indx = initial_model_indx, final_model_indx = final_model_indx)
    
    self.SBMLData = kinetics_classification.KineticAnalyzer(dataSet = "biomodels",
    model_indices=model_indices)

  def testClassification1(self):
    # Test all the column names
    if IGNORE_TEST:
      return    
    test = all(item in self.df_classification.columns for item in kinetics_classification.COLUMN_NAME_df_classification)
    self.assertTrue(test)

  def testClassification2(self):
    # Test whether there is at least one row
    if IGNORE_TEST:
      return    
    self.assertTrue(len(self.df_classification.index)>0) 

  def testClassification3(self):
    # Test all the elements of df_classification are lists of strings
    if IGNORE_TEST:
      return    
    list_classification = []
    for i in range(len(self.df_classification.columns)):
      list_classification += self.df_classification.iloc[:,i].tolist()
    test = all(isinstance(item, str) for item in list_classification)
    self.assertTrue(test)

  def testGenStat1(self):
    # Test all the column names
    if IGNORE_TEST:
      return    
    test = all(item in self.df_gen_stat.columns for item in kinetics_classification.COLUMN_NAME_df_gen_stat)
    self.assertTrue(test)

  def testGenStat2(self):
    # Test whether there is at least one row
    if IGNORE_TEST:
      return    
    self.assertTrue(len(self.df_gen_stat.index)>0) 

  def testGenStat3(self):
    # Test column 'Classification Names' of df_gen_stat is a list of string
    if IGNORE_TEST:
      return    
    list_gen_stat_classification = self.df_gen_stat[kinetics_classification.COLUMN_NAME_df_gen_stat[0]].tolist()
    test = all(isinstance(item, str) for item in list_gen_stat_classification)
    self.assertTrue(test) 

  def testGenStat4(self):
    # Test column 'Percentage' of df_gen_stat is a list of floating numbers
    if IGNORE_TEST:
      return    
    list_gen_stat_percentage = self.df_gen_stat[kinetics_classification.COLUMN_NAME_df_gen_stat[1]].tolist()
    test = all(isinstance(item, float) for item in list_gen_stat_percentage)
    self.assertTrue(test) 

  def testGenStat5(self):
    # Test column 'Percentage' of df_gen_stat does not have nan values
    if IGNORE_TEST:
      return    
    list_gen_stat_percentage = self.df_gen_stat[kinetics_classification.COLUMN_NAME_df_gen_stat[1]].tolist()
    test = any(math.isnan(item) for item in list_gen_stat_percentage)
    self.assertFalse(test)

  def testMolStat1(self):
    # Test all the column names
    if IGNORE_TEST:
      return    
    test = all(item in self.df_mol_stat.columns for item in kinetics_classification.COLUMN_NAME_df_mol_stat)
    self.assertTrue(test)
 
  def testMolStat2(self):
    # Test whether there is at least one row
    if IGNORE_TEST:
      return    
    self.assertTrue(len(self.df_mol_stat.index)>0) 

  def testMolStat3(self):
    # Test column 'SBMLid' of df_mol_stat is a list of string
    if IGNORE_TEST:
      return    
    list_gen_stat_classification = self.df_mol_stat[kinetics_classification.COLUMN_NAME_df_mol_stat[0]].tolist()
    test = all(isinstance(item, str) for item in list_gen_stat_classification)
    self.assertTrue(test)

  def testMolStat4(self):
    # Test column 'Reaction#' of df_mol_stat is a list of integers
    if IGNORE_TEST:
      return    
    list_gen_stat_classification = self.df_mol_stat[kinetics_classification.COLUMN_NAME_df_mol_stat[1]].tolist()
    test = all(isinstance(item, int) for item in list_gen_stat_classification)
    self.assertTrue(test)

  def testMolStat5(self):
    # Test columns other than 'SBMLid' and 'Reaction#' of df_mol_stat are lists of floating numbers
    if IGNORE_TEST:
      return    
    list_gen_stat_others = []
    for i in range(2,len(self.df_mol_stat.columns)):
      list_gen_stat_others += self.df_mol_stat.iloc[:,i].tolist()
    test = all(isinstance(item, float) for item in list_gen_stat_others)
    self.assertTrue(test)

  def testMolStat6(self):
    # Test whether all the numbers are not nan values
    if IGNORE_TEST:
      return    
    list_gen_stat_others = []
    for i in range(1,len(self.df_mol_stat.columns)):
      list_gen_stat_others += self.df_mol_stat.iloc[:,i].tolist()
    test = any(math.isnan(item) for item in list_gen_stat_others)
    self.assertFalse(test)

  def testMolStat7(self):
    # Test the sum of percentage of all types for each model in df_mol_stat is always no less than 1
    if IGNORE_TEST:
      return    
    list_gen_stat_others = []
    flag = 1
    for i in range(len(self.df_mol_stat)):
      sum = 0
      list_gen_stat_others += self.df_mol_stat.iloc[i,2:len(self.df_mol_stat.columns)].tolist()
      for j in range(len(list_gen_stat_others)):
        sum += list_gen_stat_others[j]
      if sum < 1:
        flag = 0
    self.assertTrue(flag)

  def testGenStatPR1(self):
    # Test all the column names
    if IGNORE_TEST:
      return    
    test = all(item in self.df_gen_stat_PR.columns for item in kinetics_classification.COLUMN_NAME_df_gen_stat[0:-2])
    self.assertTrue(test)

  def testGenStatPR2(self):
    # Test whether there is at least one row
    if IGNORE_TEST:
      return    
    self.assertTrue(len(self.df_gen_stat_PR.index)>0) 

  def testGenStatPR3(self):
    # Test column 'Percentage' of df_gen_stat is a list of floating numbers
    if IGNORE_TEST:
      return    
    list_gen_stat_percentage = self.df_gen_stat_PR[kinetics_classification.COLUMN_NAME_df_gen_stat[1]].tolist()
    test = all(isinstance(item, float) for item in list_gen_stat_percentage)
    self.assertTrue(test) 

  def testTable1(self):
    # Test all the column names
    if IGNORE_TEST:
      return    
    self.assertTrue(all(item in self.df_table_PR.columns for item in ["R = 0", "R = 1", "R = 2", "R > 2"]))
    self.assertTrue(all(item in self.df_table_PR.index for item in ["P = 0", "P = 1", "P = 2", "P > 2"]))

  def testTable2(self):
    # Test whether there is at least one row
    if IGNORE_TEST:
      return    
    self.assertTrue(len(self.df_table_PR.index)==4) 

  def testTablePerMol1(self):
    # Test all the column names
    if IGNORE_TEST:
      return    
    self.assertTrue(all(item in self.df_table_PR_per_model.columns for item in ["R = 0", "R = 1", "R = 2", "R > 2"]))
    self.assertTrue(all(item in self.df_table_PR_per_model.index for item in ["P = 0", "P = 1", "P = 2", "P > 2"]))

  def testTablePerMol2(self):
    # Test whether there is at least one row
    if IGNORE_TEST:
      return    
    self.assertTrue(len(self.df_table_PR_per_model.index)==4) 

    
  def testBiomodelNonCount1(self):
    # Test biomodel_non_count is an integer
    if IGNORE_TEST:
      return    
    test = isinstance(self.biomodel_non_count, int)
    self.assertTrue(test)

  def testGetKineticLawDistribution1(self):
    # Test getKineticLawDistribution() column names
    if IGNORE_TEST:
      return 
    df_temp = self.SBMLData.getKineticLawDistribution(fileName = "")
    test = all(item in df_temp.columns for item in ['Classifications', 'Percentage', 
    'Percentage per model', 'Percentage per model standard error'])
    self.assertTrue(test)

  def testGetKineticLawDistribution2(self):
    # Test getKineticLawDistribution() if there is at least one row
    if IGNORE_TEST:
      return 
    df_temp = self.SBMLData.getKineticLawDistribution(fileName = "")
    self.assertTrue(len(df_temp.index)>0)

  def testGetKineticLawDistribution3(self):
    # Test column 'Percentage' a list of floating numbers
    if IGNORE_TEST:
      return 
    df_temp = self.SBMLData.getKineticLawDistribution(fileName = "")
    list_percentage = df_temp['Percentage'].tolist()
    test = all(isinstance(item, float) for item in list_percentage)
    self.assertTrue(test) 

  def testGetKineticLawDistribution4(self):
    # Test column 'Percentage' does not have nan values
    if IGNORE_TEST:
      return 
    df_temp = self.SBMLData.getKineticLawDistribution(fileName = "")
    list_percentage = df_temp['Percentage'].tolist()
    test = any(math.isnan(item) for item in list_percentage)
    self.assertFalse(test)

  def testTopFrequentKineticLawDistribution(self):
    # Test TopFrequentKineticLawType()
    if IGNORE_TEST:
      return 
    self.assertTrue(self.SBMLData.TopFrequentKineticLawType() == ['ZERO', 'UNDR', 'NA'])

  def testGetKineticLawDistributionPerMassTransfer1(self):
    # Test getKineticLawDistributionPerMassTransfer() column names
    if IGNORE_TEST:
      return

    df_temp = self.SBMLData.getKineticLawDistributionPerMassTransfer(rct_num=1,prd_num=1,fileName="")
    test = all(item in df_temp.columns for item in ['Classifications', 'Percentage',
    'Percentage standard error', 'Percentage per model', 'Percentage per model standard error'])
    self.assertTrue(test)

  def testGetKineticLawDistributionPerMassTransfer2(self):
    # Test getKineticLawDistributionPerMassTransfer() if there is at least one row
    if IGNORE_TEST:
      return 
    df_temp = self.SBMLData.getKineticLawDistributionPerMassTransfer(rct_num=1,prd_num=1,fileName="")
    self.assertTrue(len(df_temp.index)>0)

  def testGetKineticLawDistributionPerMassTransfer3(self):
    # Test column 'Percentage' a list of floating numbers
    if IGNORE_TEST:
      return 
    df_temp = self.SBMLData.getKineticLawDistributionPerMassTransfer(rct_num=1,prd_num=1,fileName="")
    list_percentage = df_temp['Percentage'].tolist()
    test = all(isinstance(item, float) for item in list_percentage)
    self.assertTrue(test) 

  def testGetKineticLawDistributionPerMassTransfer4(self):
    # Test column 'Percentage' does not have nan values
    if IGNORE_TEST:
      return 
    df_temp = self.SBMLData.getKineticLawDistributionPerMassTransfer(rct_num=1,prd_num=1,fileName="")
    list_percentage = df_temp['Percentage'].tolist()
    test = any(math.isnan(item) for item in list_percentage)
    self.assertFalse(test)
  
  def testTopFrequentKineticLawDistributionPerMassTransfer(self):
    # Test TopFrequentKineticLawTypePerMassTransfer()
    if IGNORE_TEST:
      return 
    self.assertTrue(self.SBMLData.TopFrequentKineticLawTypePerMassTransfer(rct_num=1,prd_num=1) 
    == ['ZERO', 'UNDR', 'NA'])

if __name__ == '__main__':
  unittest.main()
