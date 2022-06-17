"""
Tests for kinetics_output.py
"""
# This script was written by Jin Xu and available on Github
# https://github.com/ModelEngineering/kinetics_validator
# This file includes all the tests to do the kinetics analysis.

from SBMLKinetics import kinetics_classification
from SBMLKinetics import kinetics_output
from SBMLKinetics import types
from sympy import *
import unittest
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
    = kinetics_classification._dataSetStatistics(initial_model_indx = initial_model_indx, final_model_indx = final_model_indx)
    
    self.SBMLData = kinetics_output.KineticAnalyzer(dataSet = "biomodels",
    model_indices=model_indices)

##Query Distributions
  def testGetKTypeDistribution1(self):
    # Test getKTypeDistribution() column names
    if IGNORE_TEST:
      return 
    df_temp = self.SBMLData.getKTypeDistribution()
    test = all(item in df_temp.columns for item in ['Classifications', 'Percentage', 
    'Percentage per model', 'Percentage per model standard error'])
    self.assertTrue(test)

  def testGetKTypeDistribution2(self):
    # Test getKTypeDistribution() if there is at least one row
    if IGNORE_TEST:
      return 
    df_temp = self.SBMLData.getKTypeDistribution()
    self.assertTrue(len(df_temp.index)>0)

  def testGetKTypeDistribution3(self):
    # Test column 'Percentage' a list of floating numbers
    if IGNORE_TEST:
      return 
    df_temp = self.SBMLData.getKTypeDistribution()
    list_percentage = df_temp['Percentage'].tolist()
    test = all(isinstance(item, float) for item in list_percentage)
    self.assertTrue(test) 

  def testGetKTypeDistribution4(self):
    # Test column 'Percentage' does not have nan values
    if IGNORE_TEST:
      return 
    df_temp = self.SBMLData.getKTypeDistribution()
    list_percentage = df_temp['Percentage'].tolist()
    test = any(math.isnan(item) for item in list_percentage)
    self.assertFalse(test)

  def testGetKTypeDistributionPerMType1(self):
    # Test getKTypeDistributionPerMType() column names
    if IGNORE_TEST:
      return

    df_temp = self.SBMLData.getKTypeDistributionPerMType(M_type = types.M_type(1,1))
    test = all(item in df_temp.columns for item in ['Classifications', 'Percentage',
    'Percentage standard error', 'Percentage per model', 'Percentage per model standard error'])
    self.assertTrue(test)

  def testGetKTypeDistributionPerMType2(self):
    # Test getKTypeDistributionPerMType() if there is at least one row
    if IGNORE_TEST:
      return 
    df_temp = self.SBMLData.getKTypeDistributionPerMType(M_type = types.M_type(1,1))
    self.assertTrue(len(df_temp.index)>0)

  def testGetKTypeDistributionPerMType3(self):
    # Test column 'Percentage' a list of floating numbers
    if IGNORE_TEST:
      return 
    df_temp = self.SBMLData.getKTypeDistributionPerMType(M_type = types.M_type(1,1))
    list_percentage = df_temp['Percentage'].tolist()
    test = all(isinstance(item, float) for item in list_percentage)
    self.assertTrue(test) 

  def testGetKTypeDistributionPerMType4(self):
    # Test column 'Percentage' does not have nan values
    if IGNORE_TEST:
      return 
    df_temp = self.SBMLData.getKTypeDistributionPerMType(M_type = types.M_type(1,1))
    list_percentage = df_temp['Percentage'].tolist()
    test = any(math.isnan(item) for item in list_percentage)
    self.assertFalse(test)

  def testGetMTypeDistribution1(self):
    # Test getMTypeDistribution() column names
    if IGNORE_TEST:
      return

    df_temp = self.SBMLData.getMTypeDistribution()
    test = all(item in df_temp.columns for item in ['R = 0', 'R = 1',
    'R = 2', 'R > 2'])
    self.assertTrue(test)

  def testGetMTypeDistribution2(self):
    # Test getMTypeDistribution() if there is at least one row
    if IGNORE_TEST:
      return 
    df_temp = self.SBMLData.getMTypeDistribution()
    self.assertTrue(len(df_temp.index)>0)

  def testGetMTypeDistribution3(self):
    # Test getMTypeDistribution() row names
    if IGNORE_TEST:
      return

    df_temp = self.SBMLData.getMTypeDistribution()
    test = all(item in df_temp.index for item in ['P = 0', 'P = 1',
    'P = 2', 'P > 2'])
    self.assertTrue(test)

  def testGetMTypeDistributionPerModel1(self):
    # Test getMTypeDistributionPerModel() column names
    if IGNORE_TEST:
      return

    df_temp = self.SBMLData.getMTypeDistributionPerModel()
    test = all(item in df_temp.columns for item in ['R = 0', 'R = 1',
    'R = 2', 'R > 2'])
    self.assertTrue(test)

  def testGetMTypeDistributionPerModel2(self):
    # Test getMTypeDistributionPerModel() if there is at least one row
    if IGNORE_TEST:
      return 
    df_temp = self.SBMLData.getMTypeDistributionPerModel()
    self.assertTrue(len(df_temp.index)>0)

  def testGetMTypeDistributionPerModel3(self):
    # Test getMTypeDistributionPerModel() row names
    if IGNORE_TEST:
      return

    df_temp = self.SBMLData.getMTypeDistributionPerModel()
    test = all(item in df_temp.index for item in ['P = 0', 'P = 1',
    'P = 2', 'P > 2'])
    self.assertTrue(test)

  def testGetBasicStatisticsInfo(self):
    # Test the basic statistics functions
    if IGNORE_TEST:
      return 
    self.assertTrue(self.SBMLData.getNumBiomodelsAnalyzed() == 1)
    self.assertTrue(self.SBMLData.getNumRxnsAnalyzed() == 3)

##Query Elements
  def testGetTopKType(self):
    # Test getTopKType()
    if IGNORE_TEST:
      return 
    self.assertTrue(self.SBMLData.getTopKType()[0].K_type_str == 'ZERO')
  
  def testGetKTypeProb(self):
    # Test getKTypeProb()
    if IGNORE_TEST:
      return 
    self.assertTrue(self.SBMLData.getKTypeProb(K_type = types.K_type("NA")) == 1./3)

  def testGetTopKTypePerMType(self):
    # Test getTopKTypePerMType()
    if IGNORE_TEST:
      return 
    self.assertTrue(self.SBMLData.getTopKTypePerMType(M_type = types.M_type(1,1)) 
    == ['ZERO', 'UNDR', 'NA'])

  def testGetKTypeProbPerMType(self):
    # Test getKTypeProbPerMType()
    if IGNORE_TEST:
      return 
    self.assertTrue(self.SBMLData.getKTypeProbPerMType(M_type = types.M_type(1,1), \
      K_type=types.K_type("NA")) 
    == 1./3)

  def testGetTopMType(self):
    # Test getTopMType()
    if IGNORE_TEST:
      return 

    self.assertTrue(self.SBMLData.getTopMType()[0].rct_num == 1)
  
  def testGetMTypeProb(self):
    # Test getTopMTypeProb()
    if IGNORE_TEST:
      return 
    self.assertTrue(self.SBMLData.getMTypeProb(M_type = types.M_type(1,1)) 
    == 1.)


##Presentations (not applicable)
##_functions are not applicable

if __name__ == '__main__':
  unittest.main()
