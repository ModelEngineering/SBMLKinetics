"""
Tests for Kinetic Law
"""
from SBMLKinetics.common import constants as cn
from SBMLKinetics.common import kinetic_law
from SBMLKinetics.common.simple_sbml import SimpleSBML
from SBMLKinetics.common.kinetic_law import KineticLaw
from tests.common import helpers

#from sympy import *

import copy
import libsbml
import numpy as np
import unittest


IGNORE_TEST = False
IS_PLOT = False
NUM_LAW = 3
NUM_LAW_2 = 17
NUM_LAW_3 = 7
NUM_LAW_5 = 9
NUM_LAW_43 = 7
NUM_LAW_239 = 45


class MockFunctionDefinition():
  """Used to mock function defintions."""

  def __init__(self, id, argument_names, body):
    self.id = id
    self.argument_names = argument_names
    self.body = body

  def __repr__(self):
    argument_call = ",".join(self.argument_names)
    call_str = "%s(%s)" % (self.id, argument_call)
    return "%s: %s" % (call_str, self.body)


#############################
# Tests
#############################
class TestKineticLaw(unittest.TestCase):

  def setUp(self):
    self.simple = helpers.getSimple_BIOMD6()
    self.laws = [self.simple.reactions[i].kinetic_law for i in range(NUM_LAW)]

    self.simple_2 = helpers.getSimple_BIOMD2()
    self.laws_2 = [self.simple_2.reactions[i].kinetic_law for i in range(NUM_LAW_2)]

    self.simple_3 = helpers.getSimple_BIOMD3()
    self.laws_3 = [self.simple_3.reactions[i].kinetic_law for i in range(NUM_LAW_3)]

    self.simple_5 = helpers.getSimple_BIOMD5()
    self.laws_5 = [self.simple_5.reactions[i].kinetic_law for i in range(NUM_LAW_5)]

    self.simple_43 = helpers.getSimple_BIOMD43()
    self.laws_43 = [self.simple_43.reactions[i].kinetic_law for i in range(NUM_LAW_43)]

    self.simple_239 = helpers.getSimple_BIOMD239()
    self.laws_239 = [self.simple_239.reactions[i].kinetic_law for i in range(NUM_LAW_239)]

  def testConstructor(self):
    if IGNORE_TEST:
      return
    self.assertTrue(
         isinstance(self.laws[1], KineticLaw))
    trues = [isinstance(s, str) for s in self.laws[1].symbols]
    self.assertTrue(all(trues))

  def testSymbol(self):
    if IGNORE_TEST:
      return
    def checkSubset(subset,superset):
      true = set(subset).issubset(set(superset))
      self.assertTrue(true)
    # const only
    checkSubset(['kappa'],self.laws[0].symbols)
    # simple kinetic law k*A
    checkSubset(['k6','u'],self.laws[1].symbols)
    # complex kinetic law with function pow()
    checkSubset(['z','u','k4','k4prime'],self.laws[2].symbols)

  def testExpandFormula1(self):
    # Replacement for SBML reactions
    if IGNORE_TEST:
      return
    simple = helpers.getSimple_BIOMD56()
    kinetic_law = None
    for fd in simple.function_definitions:
      for reaction in simple.reactions:
        if fd.id in reaction.kinetic_law.formula:
          kinetic_law = KineticLaw(reaction.kinetic_law.libsbml_kinetics,
              reaction, function_definitions=simple.function_definitions)
          break
      if kinetic_law is not None:
        break
    if kinetic_law is None:
      raise RuntimeError("Did not find an embedded function.")
    self.assertIsNotNone(kinetic_law.expanded_formula)
    kinetic_law_arguments = ["Vaiep", "Jaiep", "1", "IE"]
    kinetic_law.expandFormula(simple.function_definitions)
    for argument in kinetic_law_arguments:
      self.assertTrue(argument in kinetic_law.expanded_formula)

  def mkKineticLawWithFormula(self, formula):
    kinetic_law = KineticLaw(self.laws[0].libsbml_kinetics, None)
    kinetic_law.formula = formula
    return kinetic_law

  def testExpandFormula2(self):
    # Simple replacement
    if IGNORE_TEST:
      return
    function_definitions = [
        MockFunctionDefinition("aa", ["x", "y"], "x + y")
        ]
    kinetic_law = self.mkKineticLawWithFormula("aa(1, 2)")
    kinetic_law.expandFormula(function_definitions)
    self.assertEqual(kinetic_law.expanded_formula, "1 + 2")

  def testExpandFormula4(self):
    # Test replacement with other embedded expressions
    if IGNORE_TEST:
      return
    function_definitions = [
        MockFunctionDefinition("aa", ["x", "y"], "x + y")
        ]
    kinetic_law = self.mkKineticLawWithFormula("exp(4) + aa(1, 2) + sin(3)")
    kinetic_law.expandFormula(function_definitions)
    self.assertEqual(kinetic_law.expanded_formula, "exp(4) + 1 + 2 + sin(3)")

  def testExpandFormula4(self):
    # Test replacing multiple functions
    if IGNORE_TEST:
      return
    function_definitions = [
        MockFunctionDefinition("aa", ["x", "y"], "x + y"),
        MockFunctionDefinition("bb", ["x", "y"], "x*y")
        ]
    kinetic_law = self.mkKineticLawWithFormula("2 + aa(1, 2) + bb(x, z)")
    kinetic_law.expandFormula(function_definitions)
    self.assertEqual(kinetic_law.expanded_formula, "2 + 1 + 2 + x*z")

  def testExpandFormula5(self):
    # Test recursive replacements
    if IGNORE_TEST:
      return
    function_definitions = [
        MockFunctionDefinition("aa", ["x", "y"], "x + y"),
        MockFunctionDefinition("bb", ["x", "y"], "bb + aa(x, y)"),
        MockFunctionDefinition("cc", ["x", "y"], "cc + bb(x, y)")
        ]
    kinetic_law = self.mkKineticLawWithFormula("kl + cc(1, 2)")
    kinetic_law.expandFormula(function_definitions)
    self.assertEqual(kinetic_law.expanded_formula, "kl + cc + bb + 1 + 2")

  def testMkSymbolExpression(self):
    # Test replace '^' with '**'
    if IGNORE_TEST:
      return
    function_definitions = [
        MockFunctionDefinition("aa", ["x"], "x^3")
        ]
    kinetic_law = self.mkKineticLawWithFormula("aa(2) + 4^2")
    kinetic_law.mkSymbolExpression(function_definitions)
    self.assertEqual(kinetic_law.expanded_formula, "2**3 + 4**2")   

  def testIsZerothOrder1(self):
    # Test ZerothOrder success case
    if IGNORE_TEST:
      return    
    kinetic_law = self.laws[0]
    #kinetic_law: "kappa"
    species_in_kinetic_law = []
    kwargs = {"species_in_kinetic_law": species_in_kinetic_law}
    test = kinetic_law.isZerothOrder(**kwargs)
    self.assertTrue(test) 

  def testIsZerothOrder2(self):
    # Test ZerothOrder Failure case
    if IGNORE_TEST:
      return
    kinetic_law = self.laws[1]
    #kinetic_law: "k6 * u"
    species_in_kinetic_law = ['u']
    kwargs = {"species_in_kinetic_law": species_in_kinetic_law}
    test = kinetic_law.isZerothOrder(**kwargs)
    self.assertFalse(test)

  def testIsPowerTerms1(self):
    # Test Kinetics with power terms success case
    if IGNORE_TEST:
      return    
    kinetic_law = self.laws[2]
    #kinetic_law: "k4 * z * (k4prime / k4 + pow(u, 2))"
    kinetics =  "k4 * z * (k4prime / k4 + pow(u, 2))"
    try:
      kinetics_sim = str(simplify(kinetics))
    except:
      kinetics_sim = kinetics
    kwargs = {"kinetics": kinetics, "kinetics_sim": kinetics_sim}
    test = kinetic_law.isPowerTerms(**kwargs)
    self.assertTrue(test) 

  def testIsPowerTerms2(self):
    # Test Kinetics with Power terms Failure case
    if IGNORE_TEST:
      return
    kinetic_law = self.laws[1]
    #kinetic_law: "k6 * u"
    kinetics = "k6 * u"
    try:
      kinetics_sim = str(simplify(kinetics))
    except:
      kinetics_sim = kinetics
    kwargs = {"kinetics": kinetics, "kinetics_sim": kinetics_sim}
    test = kinetic_law.isPowerTerms(**kwargs)
    self.assertFalse(test)

  def testIsNoPrds1(self):
    # Test No products success case
    if IGNORE_TEST:
      return    
    kinetic_law = self.laws_3[1]
    #reaction: "C->"
    product_list = []
    kwargs = {"product_list": product_list}
    test = kinetic_law.isNoPrds(**kwargs)
    self.assertTrue(test) 

  def testIsNoPrds2(self):
    # Test No products Failure case
    if IGNORE_TEST:
      return
    kinetic_law = self.laws_3[0]
    #reaction: "->C"
    product_list = ['C']
    kwargs = {"product_list": product_list}
    test = kinetic_law.isNoPrds(**kwargs)
    self.assertFalse(test)

  def testIsSinglePrd1(self):
    # Test No products success case
    if IGNORE_TEST:
      return  
    kinetic_law = self.laws_3[0]
    #reaction: "->C"
    product_list = ['C']  
    kwargs = {"product_list": product_list}
    test = kinetic_law.isSinglePrd(**kwargs)
    self.assertTrue(test) 

  def testIsSinglePrd2(self):
    # Test No products Failure case
    if IGNORE_TEST:
      return
    kinetic_law = self.laws_3[1]
    #reaction: "C->"
    product_list = []
    kwargs = {"product_list": product_list}
    test = kinetic_law.isSinglePrd(**kwargs)
    self.assertFalse(test)

  def testIsDoublePrds1(self):
    # Test No products success case
    if IGNORE_TEST:
      return    
    kinetic_law = self.laws_5[0]
    #reaction: "M->C2 + YP"
    product_list = ['C2','YP']
    kwargs = {"product_list": product_list}
    test = kinetic_law.isDoublePrds(**kwargs)
    self.assertTrue(test) 

  def testIsDoublePrds2(self):
    # Test No products Failure case
    if IGNORE_TEST:
      return
    kinetic_law = self.laws_3[0]
    #reaction: "->C"
    product_list = ['C']
    kwargs = {"product_list": product_list}
    test = kinetic_law.isDoublePrds(**kwargs)
    self.assertFalse(test)

  def testIsMulPrds1(self):
    # Test No products success case
    if IGNORE_TEST:
      return    
    kinetic_law = self.laws_239[12]
    #reaction: "Pyr + CoA + NAD_p->CO2 + Acetyl_CoA + NADH"
    product_list = ['CO2','Acetyl_CoA','NADH']
    kwargs = {"product_list": product_list}
    test = kinetic_law.isMulPrds(**kwargs)
    self.assertTrue(test) 

  def testIsMulPrds2(self):
    # Test No products Failure case
    if IGNORE_TEST:
      return
    kinetic_law = self.laws_3[0]
    #reaction: "->C"
    product_list = ['C']
    kwargs = {"product_list": product_list}
    test = kinetic_law.isMulPrds(**kwargs)
    self.assertFalse(test)

  def testIsNoRcts1(self):
    # Test No reactants success case
    if IGNORE_TEST:
      return    
    kinetic_law = self.laws_3[0]
    #reaction: "->C"
    reactant_list = []
    kwargs = {"reactant_list": reactant_list}
    test = kinetic_law.isNoRcts(**kwargs)
    self.assertTrue(test) 

  def testIsNoRcts2(self):
    # Test No reactants Failure case
    if IGNORE_TEST:
      return
    kinetic_law = self.laws_3[1]
    #reaction: "C->"
    reactant_list = ['C']
    kwargs = {"reactant_list": reactant_list}
    test = kinetic_law.isNoRcts(**kwargs)
    self.assertFalse(test)

  def testIsSingleRct1(self):
    # Test Single reactant success case
    if IGNORE_TEST:
      return    
    kinetic_law = self.laws_3[1]
    #reaction: "C->"
    reactant_list = ['C']
    kwargs = {"reactant_list": reactant_list}
    test = kinetic_law.isSingleRct(**kwargs)
    self.assertTrue(test) 

  def testIsSingleRct2(self):
    # Test Single reactant Failure case
    if IGNORE_TEST:
      return
    kinetic_law = self.laws_3[0]
    #reaction: "->C"
    reactant_list = []
    kwargs = {"reactant_list": reactant_list}
    test = kinetic_law.isSingleRct(**kwargs)
    self.assertFalse(test)

  def testIsDoubleRcts1(self):
    # Test multiple reactants success case
    if IGNORE_TEST:
      return    
    kinetic_law = self.laws_2[0]
    #reaction: "B + L->BL"
    reactant_list = ['B', 'L']
    kwargs = {"reactant_list": reactant_list}
    test = kinetic_law.isDoubleRcts(**kwargs)
    self.assertTrue(test)

  def testIsDoubleRcts2(self):
    # Test multiple reactants Failure case
    if IGNORE_TEST:
      return
    kinetic_law = self.laws_2[3]
    #reaction: "BLL->ALL"
    reactant_list = ['BLL']
    kwargs = {"reactant_list": reactant_list}
    test = kinetic_law.isDoubleRcts(**kwargs)
    self.assertFalse(test)

  def testIsMulRcts1(self):
    # Test multiple reactants success case
    if IGNORE_TEST:
      return    
    kinetic_law = self.laws_239[12]
    #reaction: "Pyr + CoA + NAD_p->CO2 + Acetyl_CoA + NADH"
    reactant_list = ['Pyr', 'CoA', 'NAD_p']
    kwargs = {"reactant_list": reactant_list}
    test = kinetic_law.isMulRcts(**kwargs)
    self.assertTrue(test)

  def testIsMulRcts2(self):
    # Test multiple reactants Failure case
    if IGNORE_TEST:
      return
    kinetic_law = self.laws_2[3]
    #reaction: "BLL->ALL"
    reactant_list = ['BLL']
    kwargs = {"reactant_list": reactant_list}
    test = kinetic_law.isMulRcts(**kwargs)
    self.assertFalse(test)

  def testIsUNDR1(self):
    # Test Uni-directional mass reaction success case
    if IGNORE_TEST:
      return    
    kinetic_law = self.laws_43[3]
    #"Y->Z; cytosol * Kf * Y "
    reactant_list = ['Y']
    kinetics =  "cytosol * Kf * Y"
    try:
      kinetics_sim = str(simplify(kinetics))
    except:
      kinetics_sim = kinetics
    species_in_kinetic_law = ['Y']
    ids_list =  ['cytosol', 'Kf', 'Y', 'Z']
    kwargs = {"reactant_list": reactant_list, "kinetics": kinetics, "kinetics_sim": kinetics_sim, \
             "species_in_kinetic_law": species_in_kinetic_law, "ids_list": ids_list}
    test = kinetic_law.isUNDR(**kwargs)
    self.assertTrue(test) 

  def testIsUNDR2(self):
    # Test Uni-directional mass reaction Failure case
    if IGNORE_TEST:
      return
    kinetic_law = self.laws_43[0]
    #"EC->Z; cytosol * (v0 + v1 * beta)"
    reactant_list = ['EC']
    kinetics =  "cytosol * (v0 + v1 * beta) "
    try:
      kinetics_sim = str(simplify(kinetics))
    except:
      kinetics_sim = kinetics
    species_in_kinetic_law = []
    ids_list = ['cytosol', 'v0', 'v1', 'beta', 'EC', 'Z']
    kwargs = {"reactant_list": reactant_list, "kinetics": kinetics, "kinetics_sim": kinetics_sim, \
             "species_in_kinetic_law": species_in_kinetic_law, "ids_list": ids_list}
    test = kinetic_law.isUNDR(**kwargs)
    self.assertFalse(test)
  
  def testIsUNMO1(self):
    # Test Uni-term with moderator success case
    if IGNORE_TEST:
      return    
    kinetic_law = self.laws_43[5]
    #"Rho->Fraction_Inactive_Channels; cytosol * Kd * pow(Z, 4) * Rho"
    reactant_list = ['Rho']
    kinetics =  "cytosol * Kd * pow(Z, 4) * Rho"
    try:
      kinetics_sim = str(simplify(kinetics))
    except:
      kinetics_sim = kinetics
    species_in_kinetic_law = ['Z','Rho']
    ids_list = ['cytosol', 'Kd', 'Z', 'Rho', 'Fraction_Inactive_Channels']
    kwargs = {"reactant_list": reactant_list, "kinetics": kinetics, "kinetics_sim": kinetics_sim, \
             "species_in_kinetic_law": species_in_kinetic_law, "ids_list": ids_list}
    test = kinetic_law.isUNMO(**kwargs)
    self.assertTrue(test) 

  def testIsUNMO2(self):
    # Test Uni-term with moderator Failure case
    if IGNORE_TEST:
      return
    kinetic_law = self.laws_43[3]
    #"Y->Z; cytosol * Kf * Y "
    reactant_list = ['Y']
    kinetics =  "cytosol * Kf * Y "
    try:
      kinetics_sim = str(simplify(kinetics))
    except:
      kinetics_sim = kinetics
    species_in_kinetic_law = ['Y']
    ids_list = ['cytosol', 'Kf', 'Y', 'Z']
    kwargs = {"reactant_list": reactant_list, "kinetics": kinetics, "kinetics_sim": kinetics_sim, \
             "species_in_kinetic_law": species_in_kinetic_law, "ids_list": ids_list}
    test = kinetic_law.isUNMO(**kwargs)
    self.assertFalse(test)

  def testIsBIDR1(self):
    # Test Bi-directional mass reaction success case
    if IGNORE_TEST:
      return    
    kinetic_law = self.laws_2[0]
    #"B + L->BL; comp1 * (kf_0 * B * L - kr_0 * BL)"
    reactant_list = ['B','L']
    product_list = ['BL']
    kinetics =  "comp1 * (kf_0 * B * L - kr_0 * BL)"
    try:
      kinetics_sim = str(simplify(kinetics))
    except:
      kinetics_sim = kinetics
    species_in_kinetic_law = ['B','L','BL']
    ids_list = ['comp1', 'kf_0', 'B', 'L', 'kr_0', 'BL']
    kwargs = {"reactant_list": reactant_list, "product_list": product_list, \
              "kinetics": kinetics, "kinetics_sim": kinetics_sim, \
              "species_in_kinetic_law": species_in_kinetic_law, "ids_list": ids_list}
    test = kinetic_law.isBIDR(**kwargs)
    self.assertTrue(test) 

  def testIsBIDR2(self):
    # Test Bi-directional mass reaction Failure case
    if IGNORE_TEST:
      return
    kinetic_law = self.laws_43[6]
    #"Fraction_Inactive_Channels->Rho; cytosol * Kr * (1 - Rho)"
    reactant_list = ['Fraction_Inactive_Channels']
    product_list = ['Rho']
    kinetics =  "cytosol * Kr * (1 - Rho)"
    try:
      kinetics_sim = str(simplify(kinetics))
    except:
      kinetics_sim = kinetics
    species_in_kinetic_law = ['Rho']
    ids_list = ['cytosol', 'Kr', 'Rho', 'Fraction_Inactive_Channels']
    kwargs = {"reactant_list": reactant_list, "product_list": product_list, \
              "kinetics": kinetics, "kinetics_sim": kinetics_sim, \
              "species_in_kinetic_law": species_in_kinetic_law, "ids_list": ids_list}
    test = kinetic_law.isBIDR(**kwargs)
    self.assertFalse(test)

  def testIsBIMO1(self):
    # Test Bi-terms with moderator success case
    if IGNORE_TEST:
      return    
    kinetic_law = self.laws_43[6]
    #"Fraction_Inactive_Channels->Rho; cytosol * Kr * (1 - Rho)"
    reactant_list = ['Fraction_Inactive_Channels']
    product_list = ['Rho']
    kinetics =  "cytosol * Kr * (1 - Rho)"
    try:
      kinetics_sim = str(simplify(kinetics))
    except:
      kinetics_sim = kinetics
    species_in_kinetic_law = ['Rho']
    ids_list =  ['cytosol', 'Kr', 'Rho', 'Fraction_Inactive_Channels']
    kwargs = {"reactant_list": reactant_list, "product_list": product_list, \
              "kinetics": kinetics, "kinetics_sim": kinetics_sim, \
              "species_in_kinetic_law": species_in_kinetic_law, "ids_list": ids_list}
    test = kinetic_law.isBIMO(**kwargs)
    self.assertTrue(test) 

  def testIsBIMO2(self):
    # Test Bi-terms with moderator Failure case
    if IGNORE_TEST:
      return
    kinetic_law = self.laws_2[0]
    #"B + L->BL; comp1 * (kf_0 * B * L - kr_0 * BL)"
    reactant_list = ['B','L']
    product_list = ['BL']
    kinetics =  "comp1 * (kf_0 * B * L - kr_0 * BL)"
    try:
      kinetics_sim = str(simplify(kinetics))
    except:
      kinetics_sim = kinetics
    species_in_kinetic_law = ['B','L','BL']
    ids_list = ['comp1', 'kf_0', 'B', 'L', 'kr_0', 'BL']
    kwargs = {"reactant_list": reactant_list, "product_list": product_list, \
              "kinetics": kinetics, "kinetics_sim": kinetics_sim, \
              "species_in_kinetic_law": species_in_kinetic_law, "ids_list": ids_list}
    test = kinetic_law.isBIMO(**kwargs)
    self.assertFalse(test)

  def testIsMM1(self):
    # Test Michaelis-Menten Kinetics success case
    if IGNORE_TEST:
      return    
    kinetic_law = self.laws_3[4]
    #"M->; cell * M * V2 * pow(K2 + M, -1) "
    kinetics =  "cell * M * V2 * pow(K2 + M, -1)"
    try:
      kinetics_sim = str(simplify(kinetics))
    except:
      kinetics_sim = kinetics
    ids_list = ['cell', 'M', 'V2', 'K2']
    species_in_kinetic_law = ['M']
    parameters_in_kinetic_law = ['cell', 'V2', 'K2']
    reactant_list = ['M']
    kwargs = {"reactant_list": reactant_list, "kinetics": kinetics, \
              "kinetics_sim": kinetics_sim, "ids_list": ids_list, \
              "species_in_kinetic_law": species_in_kinetic_law, \
              "parameters_in_kinetic_law": parameters_in_kinetic_law}
    test = kinetic_law.isMM(**kwargs)
    self.assertTrue(test) 

  def testIsMM2(self):
    # Test Michaelis-Menten Kinetics Failure case
    if IGNORE_TEST:
      return
    kinetic_law = self.laws_3[2]
    #"C->; C * cell * vd * X * pow(C + Kd, -1)"
    kinetics =  "C * cell * vd * X * pow(C + Kd, -1)"
    try:
      kinetics_sim = str(simplify(kinetics))
    except:
      kinetics_sim = kinetics
    ids_list = ['C', 'cell', 'vd', 'X', 'Kd']
    species_in_kinetic_law = ['C', 'X']
    parameters_in_kinetic_law = ['cell', 'vd', 'Kd']
    reactant_list = ['C']
    kwargs = {"reactant_list": reactant_list, "kinetics": kinetics, \
              "kinetics_sim": kinetics_sim,"ids_list": ids_list, \
              "species_in_kinetic_law": species_in_kinetic_law, \
              "parameters_in_kinetic_law": parameters_in_kinetic_law}
    test = kinetic_law.isMM(**kwargs)
    self.assertFalse(test)

  def testIsMMcat1(self):
    # Test Michaelis-Menten Kinetics-catalyzed success case
    if IGNORE_TEST:
      return    
    kinetic_law = self.laws_3[2]
    #"C->; C * cell * vd * X * pow(C + Kd, -1)"
    kinetics =  "C * cell * vd * X * pow(C + Kd, -1)"
    try:
      kinetics_sim = str(simplify(kinetics))
    except:
      kinetics_sim = kinetics
    ids_list = ['C', 'cell', 'vd', 'X', 'Kd']
    species_in_kinetic_law = ['C', 'X']
    parameters_in_kinetic_law = ['cell', 'vd', 'Kd']
    reactant_list = ['C']
    kwargs = {"reactant_list": reactant_list, "kinetics": kinetics, \
              "kinetics_sim": kinetics_sim,"ids_list": ids_list, \
              "species_in_kinetic_law": species_in_kinetic_law, \
              "parameters_in_kinetic_law": parameters_in_kinetic_law}
    test = kinetic_law.isMMcat(**kwargs)
    self.assertTrue(test) 

  def testIsMMcat2(self):
    # Test Michaelis-Menten Kinetics-catalyzed Failure case
    if IGNORE_TEST:
      return
    kinetic_law = self.laws_3[4]
    #"M->; cell * M * V2 * pow(K2 + M, -1) "
    kinetics =  "cell * M * V2 * pow(K2 + M, -1)"
    try:
      kinetics_sim = str(simplify(kinetics))
    except:
      kinetics_sim = kinetics
    ids_list = ['cell', 'M', 'V2', 'K2']
    species_in_kinetic_law = ['M']
    parameters_in_kinetic_law = ['cell', 'V2', 'K2']
    reactant_list = ['M']
    kwargs = {"reactant_list": reactant_list, "kinetics": kinetics, \
              "kinetics_sim": kinetics_sim,"ids_list": ids_list, \
              "species_in_kinetic_law": species_in_kinetic_law, \
              "parameters_in_kinetic_law": parameters_in_kinetic_law}
    test = kinetic_law.isMMcat(**kwargs)
    self.assertFalse(test)

  def testIsHill1(self):
    # Test Hill equations success case
    if IGNORE_TEST:
      return    
    kinetic_law = self.laws_43[1]
    #"Z->Y;intravesicular * (Vm2 * pow(Z, 2) / (pow(K2, 2) + pow(Z, 2)))"
    kinetics =  "intravesicular * (Vm2 * pow(Z, 2) / (pow(K2, 2) + pow(Z, 2)))"
    try:
      kinetics_sim = str(simplify(kinetics))
    except:
      kinetics_sim = kinetics
    ids_list = ['intravesicular', 'Vm2', 'Z', 'K2', 'Y']
    species_in_kinetic_law = ['Z']
    kwargs = {"kinetics_sim": kinetics_sim,"ids_list": ids_list, \
              "species_in_kinetic_law": species_in_kinetic_law}
    test = kinetic_law.isHill(**kwargs)
    self.assertTrue(test) 

  def testIsHill2(self):
    # Test Hill equations failure case
    if IGNORE_TEST:
      return    
    kinetic_law = self.laws_43[3]
    #"Y->Z; cytosol * Kf * Y "
    kinetics =  "cytosol * Kf * Y"
    try:
      kinetics_sim = str(simplify(kinetics))
    except:
      kinetics_sim = kinetics
    species_in_kinetic_law = ['Y']
    ids_list =  ['cytosol', 'Kf', 'Y', 'Z']
    kwargs = {"kinetics_sim": kinetics_sim,"ids_list": ids_list, \
              "species_in_kinetic_law": species_in_kinetic_law}
    test = kinetic_law.isHill(**kwargs)
    self.assertFalse(test) 

  def testIsFraction1(self):
    # Test Fraction success case
    if IGNORE_TEST:
      return    
    kinetic_law = self.laws_3[2]
    #kinetic_law: "C * cell * vd * X * pow(C + Kd, -1)"
    kinetics_sim = "C*X*cell*vd*pow(C + Kd, -1)"
    ids_list = ['C', 'cell', 'vd', 'X', 'Kd']
    species_in_kinetic_law = ['C', 'X'] 

    kwargs = {"kinetics_sim": kinetics_sim, "ids_list": ids_list, "species_in_kinetic_law": species_in_kinetic_law}
    test = kinetic_law.isFraction(**kwargs)
    self.assertTrue(test) 

  def testIsFraction2(self):
    # Test Fraction Failure case
    if IGNORE_TEST:
      return
    kinetic_law = self.laws_3[0]
    #kinetic_law: "cell * vi"
    kinetics_sim = "cell*vi"
    ids_list = ['cell', 'vi', 'C']
    species_in_kinetic_law = [] 

    kwargs = {"kinetics_sim": kinetics_sim, "ids_list": ids_list, "species_in_kinetic_law": species_in_kinetic_law}
    test = kinetic_law.isFraction(**kwargs)
    self.assertFalse(test)

  def testIsPolynomial1(self):
    # Test Polynomial success case
    if IGNORE_TEST:
      return    
    kinetic_law = self.laws_3[1]
    #kinetic_law: "C * cell * kd"
    kinetics_sim = "C * cell * kd"
    ids_list = ['C', 'cell', 'kd']
    species_in_kinetic_law = ['C'] 

    kwargs = {"kinetics_sim": kinetics_sim, "ids_list": ids_list, "species_in_kinetic_law": species_in_kinetic_law}
    test = kinetic_law.isPolynomial(**kwargs)
    self.assertTrue(test) 

  def testIsPolynomial2(self):
    # Test Polynomial Failure case
    if IGNORE_TEST:
      return
    kinetic_law = self.laws_3[0]
    #kinetic_law: "cell * vi"
    kinetics_sim = "cell*vi"
    ids_list = ['cell', 'vi', 'C']
    species_in_kinetic_law = [] 

    kwargs = {"kinetics_sim": kinetics_sim, "ids_list": ids_list, "species_in_kinetic_law": species_in_kinetic_law}
    test = kinetic_law.isPolynomial(**kwargs)
    self.assertFalse(test)

if __name__ == '__main__':
  unittest.main()

