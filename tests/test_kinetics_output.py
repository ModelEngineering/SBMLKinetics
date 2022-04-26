"""
Tests for kinetics_classification.py
"""
from SBMLKinetics import kinetics_classification
from SBMLKinetics import kinetics_output
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
