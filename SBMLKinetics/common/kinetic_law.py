'''Provides Information on SBML Kinetics Laws'''
# This script was written by Jin Xu and available on Github
# https://github.com/SunnyXu/SBMLKinetics
# This file includes all the functions of types to classify the kinetics.


from SBMLKinetics.common import constants as cn
from SBMLKinetics.common import util
from SBMLKinetics.common import exceptions
from SBMLKinetics.common import msgs
import sympy
from sympy import symbols
from sympy import core

import collections #use set to compare two lists
import numpy as np
import re # Extract substrings between brackets

MAX_RECURSION = 5 # Maximum number for iteration function expansions

class KineticLaw(object):

  def __init__(self, libsbml_kinetics, reaction, function_definitions=None):
    """
    :param libsbml.KineticLaw libsbml_kinetics:
    :param function_definitions list-FunctionDefinition:
    """
    # libsbml object for kinetics
    self.libsbml_kinetics = libsbml_kinetics
    # String version of chemical formula
    try:
      self.formula = self.libsbml_kinetics.getFormula()
    except:
      self.formula = ""
    # Reaction for the kinetics law
    self.reaction = reaction
    # Parameters and chemical species
    # self.symbols = self._getSymbols()
    try:
      self.symbols = self._getSymbols()
    except Exception:
      self.symbols = []
    # Expanded kinetic formula (remove embedded functions)
    if function_definitions is None:
      self.expanded_formula = None
    else:
      self.expandFormula(function_definitions)
    self.expression_formula = None  # valid symPy expression string

  def __repr__(self):
    return self.formula

  def expandFormula(self, function_definitions):
    """
    Expands the kinetics formula, replacing function definitions
    with their body.

    Parameters
    ----------
    function_definitions: list-FunctionDefinition
    """
    self.expanded_formula = self._expandFormula(self.formula, function_definitions)

  def mkSymbolExpression(self, function_definitions):
    """
    Creates a string that can be processed by sympy.
    
    Parameters
    -------
    function_definitions: list-FunctionDefinition
    
    Returns
    -------
    str
    """
    #if self.expanded_formula is/is not?? None:
    self.expandFormula(function_definitions)
    self.expression_formula = str(self.expanded_formula)
    self.expanded_formula = self.expression_formula.replace("^","**")
    return self.expanded_formula

  def isZerothOrder(self, **kwargs):
    """
    Check whether the reaction with a kinetic law belongs to the type of Zeroth Order 
    Zeroth order classification rule: if there are no species in the kinetics
    
    Parameters
    -------
    species_in_kinetic_law: list-species in the kinetics
    
    Returns
    -------
    True or False
    """
    species_in_kinetic_law = kwargs["species_in_kinetic_law"]
    return self._numSpeciesInKinetics(species_in_kinetic_law) == 0  
  
  def isPowerTerms(self, **kwargs):
    """
    Check whether the reaction with a kinetic law belongs to the type of Kinetics with power terms
    Kinetics with power terms classification rule: if there is pow() or ** inside the kinetics, 
    except the pow(,-1) case as the possible Michaelis–Menten kinetics

    Parameters
    -------
    **kwargs: dictionary-keyword arguments
    kinetics: string-kinetics
    kinetics_sim: string-simplified kinetics
    
    Returns
    -------
    True or False
    """
    kinetics = kwargs["kinetics"]
    kinetics_sim = kwargs["kinetics_sim"]
    return self._powerInKinetics(kinetics, kinetics_sim) 

  def isNoPrds(self, **kwargs):
    """
    Tests whether the reaction belongs to the type of No Products
    No products classification rule: if there are no products

    Parameters
    -------
    **kwargs: dictionary-keyword arguments
    product_list: list-products of the reaction
    
    Returns
    -------
    True or False
    """
    product_list = kwargs["product_list"]
    return self._numOfPrds(product_list) == 0
  
  def isSinglePrd(self, **kwargs):
    """
    Tests whether the reaction belongs to the type of Single Product
    Single product classification rule: if there is a single product

    Parameters
    -------
    **kwargs: dictionary-keyword arguments
    product_list: list-products of the reaction
    
    Returns
    -------
    True or False
    """
    product_list = kwargs["product_list"]
    return self._numOfPrds(product_list) == 1

  def isDoublePrds(self, **kwargs):
    """
    Tests whether the reaction belongs to the type of Double Products
    Double products classification rule: if there are double products

    Parameters
    -------
    **kwargs: dictionary-keyword arguments
    product_list: list-products of the reaction
    
    Returns
    -------
    True or False
    """
    product_list = kwargs["product_list"]
    return self._numOfPrds(product_list) == 2

  def isMulPrds(self, **kwargs):
    """
    Tests whether the reaction belongs to the type of Multiple Products
    Multiple products classification rule: if there is are more than two products

    Parameters
    -------
    **kwargs: dictionary-keyword arguments
    product_list: list-products of the reaction
    
    Returns
    -------
    True or False
    """
    product_list = kwargs["product_list"]
    return self._numOfPrds(product_list) > 2

  def isNoRcts(self, **kwargs):
    """
    Tests whether the reaction belongs to the type of No Reactants
    No reactants classifcation rule: if there are no reactants

    Parameters
    -------
    **kwargs: dictionary-keyword arguments
    reactant_list: list-reactants of the reaction
    
    Returns
    -------
    True or False
    """
    reactant_list = kwargs["reactant_list"]
    return self._numOfRcts(reactant_list) == 0

  def isSingleRct(self, **kwargs):
    """
    Tests whether the reaction belongs to the type of Single Reactant
    Single reactant classification rule: if there is only one reactant
    
    Parameters
    -------
    **kwargs: dictionary-keyword arguments
    reactant_list: list-reactants of the reaction
    
    Returns
    -------
    True or False
    """
    reactant_list = kwargs["reactant_list"]
    return self._numOfRcts(reactant_list) == 1

  def isDoubleRcts(self, **kwargs):
    """
    Tests whether the reaction belongs to the type of Double Reactants
    Double reactants classification rule: if there is two reactants
    
    Parameters
    -------
    **kwargs: dictionary-keyword arguments
    reactant_list: list-reactants of the reaction
    
    Returns
    -------
    True or False
    """
    reactant_list = kwargs["reactant_list"]
    return self._numOfRcts(reactant_list) == 2

  def isMulRcts(self, **kwargs):
    """
    Tests whether the reaction belongs to the type of Multiple Reactants
    Multiple reactants classification rule: if there are more than two reactants
    
    Parameters
    -------
    **kwargs: dictionary-keyword arguments
    reactant_list: list-reactants of the reaction
    
    Returns
    -------
    True or False
    """
    reactant_list = kwargs["reactant_list"] 
    return self._numOfRcts(reactant_list) > 2


  def isUNDR(self, **kwargs):
    """
    Tests whether the reaction belongs to the type of uni-directional mass reaction
    Uni-directional mass reaction classification rule: 
    1) Kinetics is a single product of terms or with an additional const term
    2) The species inside the kinetics are only reactants
    
    Parameters
    -------
    **kwargs: dictionary-keyword arguments
    reactant_list: list-reactants of the reaction
    kinetics: string-kinetics
    kinetics_sim: string-simplified kinetics
    species_in_kinetic_law: list-species in the kinetics
    ids_list: list-id list including all the ids in kinetics, reactants and products
    
    Returns
    -------
    True or False
    """
    reactant_list = kwargs["reactant_list"]
    kinetics = kwargs["kinetics"]
    kinetics_sim = kwargs["kinetics_sim"]
    species_in_kinetic_law = kwargs["species_in_kinetic_law"]
    ids_list = kwargs["ids_list"]

    flag = False
    if self._isSingleProductOfTerms(kinetics, kinetics_sim) == True \
      and self._SpecsInKineticsAllRcts(species_in_kinetic_law, reactant_list) == True:
      flag = True
    if len(species_in_kinetic_law) == 1 and species_in_kinetic_law == reactant_list:
      if kinetics.count(species_in_kinetic_law[0]) == 1:
        flag = True
      elif kinetics_sim.count(species_in_kinetic_law[0]) == 1:
        flag = True 
    try:
      eq = self._numeratorDenominator(kinetics_sim, ids_list)
      if len(species_in_kinetic_law) > 0:
        for i in range(len(species_in_kinetic_law)):
          if species_in_kinetic_law[i] in eq[1]:
            flag = False
    except:
      pass

    return flag

  def isUNMO(self, **kwargs):
    """
    Tests whether the reaction belongs to the type of uni-term with moderator
    Uni-term with moderator classification rule: 
    1) Kinetics is a single product of terms
    2) The species inside the kinetics are not only reactants
    
    Parameters
    -------
    **kwargs: dictionary-keyword arguments
    reactant_list: list-reactants of the reaction
    kinetics: string-kinetics
    kinetics_sim: string-simplified kinetics
    species_in_kinetic_law: list-species in the kinetics
    ids_list: list-id list including all the ids in kinetics, reactants and products
    
    Returns
    -------
    True or False
    """
    reactant_list = kwargs["reactant_list"]
    kinetics = kwargs["kinetics"]
    kinetics_sim = kwargs["kinetics_sim"]
    species_in_kinetic_law = kwargs["species_in_kinetic_law"]
    ids_list = kwargs["ids_list"]

    flag = False
    if self._isSingleProductOfTerms(kinetics, kinetics_sim) == True \
      and self._SpecsInKineticsAllRcts(species_in_kinetic_law, reactant_list) == False\
      and self._numSpeciesInKinetics(species_in_kinetic_law) != 0:
      flag = True
    if len(species_in_kinetic_law) == 1 and species_in_kinetic_law != reactant_list\
      and self._isDiffOfTwoProductsOfTerms(kinetics, kinetics_sim) == False:
      if kinetics.count(species_in_kinetic_law[0]) == 1:
        flag = True
      elif kinetics_sim.count(species_in_kinetic_law[0]) == 1:
        flag = True 
    try:
      eq = self._numeratorDenominator(kinetics_sim, ids_list)
      if len(species_in_kinetic_law) > 0:
        for i in range(len(species_in_kinetic_law)):
          if species_in_kinetic_law[i] in eq[1]:
            flag = False
    except:
      pass
    
    return flag

  def isBIDR(self, **kwargs):
    """
    Tests whether the reaction belongs to the type of bi-directional mass reaction
    Bi-directional mass reactionclassification rule: 
    1) Kinetics is the difference of two product of terms, with species in each term
    2) The first term before - includes all the reactants
       while the second term after - includes all the products
    
    Parameters
    -------
    **kwargs: dictionary-keyword arguments
    reactant_list: list-reactants of the reaction
    product_list: list-products of the reaction
    kinetics: string-kinetics
    kinetics_sim: string-simplified kinetics
    species_in_kinetic_law: list-species in the kinetics
    ids_list: list-id list including all the ids in kinetics, reactants and products

    Returns
    -------
    True or False
    """

    reactant_list = kwargs["reactant_list"]
    product_list = kwargs["product_list"]
    kinetics = kwargs["kinetics"]
    kinetics_sim = kwargs["kinetics_sim"]
    species_in_kinetic_law = kwargs["species_in_kinetic_law"]
    ids_list = kwargs["ids_list"]

    flag = True
    if self._isDiffOfTwoProductsOfTerms(kinetics, kinetics_sim) == False:
      flag = False
    if self._ProductOfTermsWithAllRctsOrPrds(kinetics, kinetics_sim, species_in_kinetic_law, reactant_list, product_list) == False:
      flag = False
    try:
      eq = self._numeratorDenominator(kinetics_sim, ids_list)
      if len(species_in_kinetic_law) > 0:
        for i in range(len(species_in_kinetic_law)):
          if species_in_kinetic_law[i] in eq[1]:
            flag = False
    except:
      pass
    
    return flag

  def isBIMO(self, **kwargs):
    """
    Tests whether the reaction belongs to the type of bi-terms with moderator
    Bi-terms with moderator classification rule: 
    1) Kinetics is the difference of two product of terms
    2) The first term before - does not include all the reactants
       while the second term after - does not include all the products
    
    Parameters
    -------
    **kwargs: dictionary-keyword arguments
    reactant_list: list-reactants of the reaction
    product_list: list-products of the reaction
    kinetics: string-kinetics
    kinetics_sim: string-simplified kinetics
    species_in_kinetic_law: list-species in the kinetics
    ids_list: list-id list including all the ids in kinetics, reactants and products
    
    Returns
    -------
    True or False
    """
    reactant_list = kwargs["reactant_list"]
    product_list = kwargs["product_list"]
    kinetics = kwargs["kinetics"]
    kinetics_sim = kwargs["kinetics_sim"]
    species_in_kinetic_law = kwargs["species_in_kinetic_law"]
    ids_list = kwargs["ids_list"]

    flag = True
    if self._numSpeciesInKinetics(species_in_kinetic_law) == 0: #exclude the case of ZERO
      flag = False
    if self._isDiffOfTwoProductsOfTerms(kinetics, kinetics_sim) == False:
      flag = False
    if self._ProductOfTermsWithAllRctsOrPrds(kinetics, kinetics_sim, species_in_kinetic_law, reactant_list, product_list) == True:
      flag = False
    if len(species_in_kinetic_law) == 1 and species_in_kinetic_law == reactant_list: 
    #exclude the case of UNDR
      if kinetics.count(species_in_kinetic_law[0]) == 1:
        flag = False
      elif kinetics_sim.count(species_in_kinetic_law[0]) == 1:
        flag = False 
    try:
      eq = self._numeratorDenominator(kinetics_sim, ids_list)
      if len(species_in_kinetic_law) > 0:
        for i in range(len(species_in_kinetic_law)):
          if species_in_kinetic_law[i] in eq[1]:
            flag = False
    except:
      pass

    return flag

  def isMM(self, **kwargs):

    """
    Tests whether the reaction belongs to the type of Michaelis-Menten Kinetics
    Michaelis–Menten kinetics(inreversible) classification rule:
    assuming there are one/two/three parameters in the numerator, using "simplify" equals to
    
    Parameters
    ---- 
    **kwargs: dictionary-keyword arguments 
    kinetics: string-kinetics
    kinetics_sim: string-simplified kinetics
    ids_list: list-id list including all the ids in kinetics, reactants and products
    species_in_kinetic_law: list-species in the kinetics
    parameters_in_kinetic_law: list-parameters in the kinetics  
    reactant_list: list-reactants of the reaction
    
    Returns
    -------
    True or False
    """
  
    kinetics = kwargs["kinetics"]
    kinetics_sim = kwargs["kinetics_sim"]
    ids_list = kwargs["ids_list"]
    species_in_kinetic_law = kwargs["species_in_kinetic_law"]
    parameters_in_kinetic_law = kwargs["parameters_in_kinetic_law"]
    reactant_list = kwargs["reactant_list"]


    eq = self._numeratorDenominator(kinetics_sim, ids_list)
    flag_fr = False
    if len(species_in_kinetic_law) > 0:
      for i in range(len(species_in_kinetic_law)):
        if species_in_kinetic_law[i] in eq[1]:
          flag_fr = True

    flag = False
    if flag_fr:
      if self._numSpeciesInKinetics(species_in_kinetic_law) == 1 and self._numOfRcts(reactant_list) == 1:
        if self._MMSingleSpecInNumerator(kinetics, ids_list, parameters_in_kinetic_law, reactant_list) == True:
          flag = True
    else:
      flag = False

    return flag


  def isMMcat(self, **kwargs):
    """
    Tests whether the reaction belongs to the type of Michaelis-Menten Kinetics-catalyzed
    Michaelis–Menten kinetics(catalyzed) classification rule:
    Assuming there are no/one/two parameters in the numerator, using "simplify" equals to
    
    Parameters
    ----  
    **kwargs: dictionary-keyword arguments  
    kinetics: string-kinetics
    ids_list: list-id list including all the ids in kinetics, reactants and products
    species_in_kinetic_law: list-species in the kinetics
    parameters_in_kinetic_law: list-parameters in the kinetics
    reactant_list: list-reactants of the reaction

    Returns
    -------
    True or False
    """
      
    kinetics = kwargs["kinetics"]
    kinetics_sim = kwargs["kinetics_sim"]
    ids_list = kwargs["ids_list"]
    species_in_kinetic_law = kwargs["species_in_kinetic_law"]
    parameters_in_kinetic_law = kwargs["parameters_in_kinetic_law"]
    reactant_list = kwargs["reactant_list"]

    eq = self._numeratorDenominator(kinetics_sim, ids_list)
    flag_fr = False
    if len(species_in_kinetic_law) > 0:
      for i in range(len(species_in_kinetic_law)):
        if species_in_kinetic_law[i] in eq[1]:
          flag_fr = True

    flag = False
    if flag_fr:
      if self._numSpeciesInKinetics(species_in_kinetic_law) == 2 and self._numOfRcts(reactant_list) == 1:
        if self._MMTwoSpecInNumerator(kinetics, ids_list, parameters_in_kinetic_law, species_in_kinetic_law, reactant_list) == True:  
          flag = True
    else:
      flag = False
      
    return flag

  def isHill(self, **kwargs):

    """
    Tests whether the reaction belongs to the type of Hill equations.
    Hill equations classification rule:
    1) kinetics has to be in the format of fraction;
    2) there is only one species in the kinetics, but it does not have to be a reactant;
    3) kinetics is in the Hill format;

    Parameters
    ---- 
    **kwargs: dictionary-keyword arguments 
    kinetics_sim: string-simplified kinetics
    ids_list: list-id list including all the ids in kinetics, reactants and products
    species_in_kinetic_law: list-species in the kinetics
    
    Returns
    -------
    True or False
    """

    kinetics_sim = kwargs["kinetics_sim"]
    ids_list = kwargs["ids_list"]
    species_in_kinetic_law = kwargs["species_in_kinetic_law"]

    eq = self._numeratorDenominator(kinetics_sim, ids_list)
    flag_fr = False
    if len(species_in_kinetic_law) > 0:
      for i in range(len(species_in_kinetic_law)):
        if species_in_kinetic_law[i] in eq[1]:
          flag_fr = True

    flag = False
    if flag_fr:
      if self._numSpeciesInKinetics(species_in_kinetic_law) == 1:
        if self._HillFormat(kinetics_sim, ids_list, species_in_kinetic_law) == True:
          flag = True
    else:
      flag = False


    return flag

  def isFraction(self, **kwargs):
    """
    Tests whether the reaction belongs to the type of Fraction function
    Fraction function classification rule:
    Check whether it is in the form of fraction and has species in the denominator
    
    Parameters
    ----  
    **kwargs: dictionary-keyword arguments  
    kinetics_sim: string-simplified kinetics
    ids_list: list-id list including all the ids in kinetics, reactants and products
    species_in_kinetic_law: list-species in the kinetics

    Returns
    -------
    True or False
    """
    kinetics_sim = kwargs["kinetics_sim"]
    ids_list = kwargs["ids_list"]
    species_in_kinetic_law = kwargs["species_in_kinetic_law"]

    eq = self._numeratorDenominator(kinetics_sim, ids_list)
    flag = False
    if len(species_in_kinetic_law) > 0:
      for i in range(len(species_in_kinetic_law)):
        if species_in_kinetic_law[i] in eq[1]:
          flag = True
    return flag

  def isPolynomial(self, **kwargs):
    """
    Tests whether the reaction belongs to the type of Polynomial function
    Polynomial function classification rule:
    Check whether it is in the form of polynomial
    
    Parameters
    ----  
    **kwargs: dictionary-keyword arguments  
    kinetics_sim: string-simplified kinetics
    ids_list: list-id list including all the ids in kinetics, reactants and products
    species_in_kinetic_law: list-species in the kinetics

    Returns
    -------
    True or False
    """
    kinetics_sim = kwargs["kinetics_sim"]
    ids_list = kwargs["ids_list"]
    species_in_kinetic_law = kwargs["species_in_kinetic_law"]

    flag = False
    if self._isPolynomial(kinetics_sim, ids_list) == True and len(species_in_kinetic_law) > 0:
      for i in range(len(species_in_kinetic_law)):
        if species_in_kinetic_law[i] in kinetics_sim:
          flag = True
    return flag
    
  def _numSpeciesInKinetics(self, species_in_kinetic_law):
    """
    Tests whether there is no species in the kinetic law
    
    Parameters
    -------
    species_in_kinetic_law: list-species in the kinetics

    Returns
    -------
    Integer
    """
    return len(species_in_kinetic_law)

  def _powerInKinetics(self, kinetics, kinetics_sim):
    """
    Tests whether there is power term in the kinetic law: **, "pow", but not "pow(,-1)"
    
    Parameters
    ----    
    kinetics: string-kineticse
    kinetics_sim: string-simplified kinetics
    
    Returns
    -------
    True or False
    """
    if "pow(" in kinetics and "-1)" not in kinetics:
      return True
    elif "**" in kinetics: 
      return True
    elif "**" in kinetics_sim:
      return True
    else:
      False
  
  def _numOfPrds(self, product_list):
    """
    Tests for the number of prds in the reaction
    
    Parameters
    ----    
    product_list: list-products of the reaction
    
    Returns
    -------
    Integer
    """
    return len(product_list)

  def _numOfRcts(self, reactant_list):
    """
    Tests for the number of rcts in the reaction
    
    Parameters
    ----    
    reactant_list: list-reactants of the reaction
    
    Returns
    -------
    Integer
    """
    return len(reactant_list)

  def _isSingleProductOfTerms(self, kinetics, kinetics_sim):
    """
    Tests whether the kinetics is a single product of terms
    
    Parameters
    ----    
    kinetics: string-kinetics
    kinetics_sim: string-simplified kinetics
    
    Returns
    -------
    True or False
    """
    flag = True
    if "+" in kinetics or "-" in kinetics:
      flag = False
      if "e-" in kinetics or "exp(-" in kinetics:
        flag = True
    elif "+" in kinetics_sim or "-" in kinetics_sim:
      flag = False
      if "e-" in kinetics or "exp(-" in kinetics_sim:
        flag = True
    return flag

  def _SpecsInKineticsAllRcts(self, species_in_kinetic_law, reactant_list):
    """
    Tests whether all species in kinetics are reactants
    
    Parameters
    ----    
    species_in_kinetic_law: list-species in the kinetics
    reactant_list: list-reactants of the reaction
    
    Returns
    -------
    True or False
    """
    if len(reactant_list) > 0 and collections.Counter(species_in_kinetic_law) == collections.Counter(reactant_list):
      return True
    else:
      return False

  def _isDiffOfTwoProductsOfTerms(self, kinetics, kinetics_sim):
    """
    Tests whether the kinetics is the difference between two product of terms
    
    Parameters
    ----    
    kinetics: string-kinetics
    kinetics_sim: string-simplified kinetics
    
    Returns
    -------
    True or False
    """

    flag = False
    if self._isSingleProductOfTerms(kinetics, kinetics_sim) == False:
      terms = kinetics.split("-")
      if len(terms) == 2:
        flag = True
    return flag

  def _ProductOfTermsWithAllRctsOrPrds(self, kinetics, kinetics_sim, species_in_kinetic_law, reactant_list, product_list):
    """
    Tests whether the kinetics with one/the other product terms with all reactants/products
    
    Parameters
    ----    
    kinetics: string-kinetics
    kinetics_sim: string-simplified kinetics
    species_in_kinetic_law: list-species in the kinetics
    reactant_list: list-reactants of the reaction
    product_list: list-products of the reaction

    Returns
    -------
    True or False
    """

    flag_kinetics = 1
    flag_kinetics_sim = 1
    terms = kinetics.split("-") 
    if len(terms) == 2 and len(reactant_list) > 0 and len(product_list) > 0:
      term1 = terms[0]
      term2 = terms[1]
      if collections.Counter(species_in_kinetic_law) == collections.Counter(reactant_list+product_list):
        if all(ele in term1 for ele in reactant_list) and all(ele in term2 for ele in product_list):
          flag_kinetics = 0

    terms = kinetics_sim.split("-") 
    if len(terms) == 2 and len(reactant_list) > 0 and len(product_list) > 0:
      term1 = terms[0]
      term2 = terms[1]
      if collections.Counter(species_in_kinetic_law) == collections.Counter(reactant_list+product_list):
        if all(ele in term1 for ele in reactant_list) and all(ele in term2 for ele in product_list):
          flag_kinetics_sim = 0

    if flag_kinetics*flag_kinetics_sim == 0:
      return True
    else:
      return False


  def _MMSingleSpecInNumerator(self, kinetics, ids_list, parameters_in_kinetic_law, reactant_list):
    """
    Tests whether kinetics is in the MM functional form with a single species in the numerator
    
    Parameters
    ----    
    kinetics: string-kinetics
    ids_list: list-id list including all the ids in kinetics, reactants and products
    parameters_in_kinetic_law: list-parameters in the kinetics  
    reactant_list: list-reactants of the reaction
    
    Returns
    -------
    True or False
    """

    strange_func = 0 #check if there are strang functions (i.e. delay) in kinetics
    flag = 0
    pre_symbols = ''
    for i in range(len(ids_list)):
      pre_symbols += ids_list[i]
      pre_symbols += ' '
    pre_symbols = pre_symbols[:-1] #remove the space at the end
    pre_symbols_comma = pre_symbols.replace(" ",",")
    stmt = "%s = symbols('%s')"%(pre_symbols_comma,pre_symbols)
    try: #sometimes there is "invalid syntax error"
      exec(stmt,globals())
    except: 
      strange_func = 1

    try: #check if there is strange func (i.e. delay) in kinetic law
      expr_stat = "expr = " + kinetics
      exec(expr_stat,globals())
    except:
      strange_func = 1

    if strange_func == 0:
      if (len(parameters_in_kinetic_law) >= 2):                    
        for j in range(len(parameters_in_kinetic_law)):
          for k in range(len(parameters_in_kinetic_law)):
            for l in range(len(parameters_in_kinetic_law)):
              for m in range(len(parameters_in_kinetic_law)):
                # assuming there is one parameter in the numerator
                if k != j:
                  pre_n = reactant_list[0]
                  pre_d = ' ( '
                  pre_d += reactant_list[0]
                  pre_n += ' * '
                  pre_d += ' + '
                  pre_n += parameters_in_kinetic_law[j]
                  pre_d += parameters_in_kinetic_law[k]
                  pre_d += ' ) '
                  pre = pre_n
                  pre += ' / '
                  pre += pre_d          
                  expr1_stat = "expr1 =" + pre
                  exec(expr1_stat,globals())
                  if sympy.simplify(expr1) == sympy.simplify(expr):
                    flag = 1
                    break
                  # assuming there are two parameters in the numerator
                  if len(parameters_in_kinetic_law) >= 3:
                    if l != j and l != k:
                      pre_n += ' * '
                      pre_n += parameters_in_kinetic_law[l]
                      pre = pre_n
                      pre += ' / '
                      pre += pre_d           
                      expr1_stat = "expr1 =" + pre
                      exec(expr1_stat,globals()) 
                      #exec() does not work in python function?
                      if sympy.simplify(expr1) == sympy.simplify(expr):
                        flag = 1
                        break
                      # assuming there are three parameters in the numerator
                      if len(parameters_in_kinetic_law) >= 4:
                        if m != j and m != k and m != l:
                          pre_n += ' * '
                          pre_n += parameters_in_kinetic_law[m]
                          pre = pre_n
                          pre += ' / '
                          pre += pre_d       
                          expr1_stat = "expr1 =" + pre
                          exec(expr1_stat,globals())
                          if sympy.simplify(expr1) == sympy.simplify(expr):
                            flag = 1
                            break
              if flag == 1:
                break
            if flag == 1:
              break
          if flag == 1:
            break
    
    if flag == 1:
      return True
    else: 
      return False

  def _MMTwoSpecInNumerator(self, kinetics, ids_list, parameters_in_kinetic_law, species_in_kinetic_law, reactant_list):
    """
    Tests whether kinetics is in the MM functional with a reactant and 2nd species as a product in the numerator
    
    Parameters
    ----    
    kinetics: string-kinetics
    ids_list: list-id list including all the ids in kinetics, reactants and products
    parameters_in_kinetic_law: list-parameters in the kinetics 
    species_in_kinetic_law: list-species in the kinetics
    reactant_list: list-reactants of the reaction
    
    Returns
    -------
    True or False
    """

    strange_func = 0 #check if there are strang functions (i.e. delay) in kinetics
    flag = 0
    pre_symbols = ''
    for i in range(len(ids_list)):
      pre_symbols += ids_list[i]
      pre_symbols += ' '
    pre_symbols = pre_symbols[:-1] #remove the space at the end
    pre_symbols_comma = pre_symbols.replace(" ",",")
    stmt = "%s = symbols('%s')"%(pre_symbols_comma,pre_symbols)
    try: #sometimes there is "invalid syntax error"
      exec(stmt,globals())
    except: 
      strange_func = 1

    try: #check if there is strange func (i.e. delay) in kinetic law
      expr_stat = "expr = " + kinetics
      exec(expr_stat,globals())
    except:
      strange_func = 1

    if strange_func == 0:
      if (len(parameters_in_kinetic_law) != 0):                    
        for j in range(len(parameters_in_kinetic_law)):
          for k in range(len(parameters_in_kinetic_law)):
            for l in range(len(parameters_in_kinetic_law)):
              #no parameter in the numerator
              pre_n = reactant_list[0]
              cat = [item for item in species_in_kinetic_law if item not in reactant_list][0]
              pre_n += ' * '
              pre_n += cat
              pre_d = ' ( '
              pre_d += reactant_list[0]
              pre_d += ' + '
              pre_d += parameters_in_kinetic_law[k]
              pre_d += ' ) '
              pre = pre_n
              pre += ' / '
              pre += pre_d           
              expr1_stat = "expr1 =" + pre
              exec(expr1_stat,globals())
              if sympy.simplify(expr1) == sympy.simplify(expr):
                flag = 1
                break
              # assuming there is one parameter in the numerator
              if len(parameters_in_kinetic_law) >= 2:
                if k != j:
                  pre_n += ' * '
                  pre_n += parameters_in_kinetic_law[j]
                  pre = pre_n
                  pre += ' / '
                  pre += pre_d           
                  expr1_stat = "expr1 =" + pre
                  exec(expr1_stat,globals())
                  if sympy.simplify(expr1) == sympy.simplify(expr):
                    flag = 1
                    break
                  # assuming there are two parameters in the numerator
                  if len(parameters_in_kinetic_law) >= 3:
                    if l != j and l != k:
                      pre_n += ' * '
                      pre_n += parameters_in_kinetic_law[l]
                      pre = pre_n
                      pre += ' / '
                      pre += pre_d
                      expr1_stat = "expr1 =" + pre
                      exec(expr1_stat,globals())
                      if sympy.simplify(expr1) == sympy.simplify(expr):
                        flag = 1
                        break
            if flag == 1:
              break
          if flag == 1:
            break        
    
    if flag == 1:
      return True
    else: 
      return False

  def _HillFormat(self, kinetics_sim, ids_list, species_in_kinetic_law):
    """
    Tests whether the kinetics is in the format of Hill equations.
    1) the numerator is one product of terms with the species to a power;
    2) the denomimator is the sum of two product of terms, one of which does not include the species
       and the other one of which include the species to the same power as the numerator.
    
    Parameters
    ----    
    kinetics_sim: string-simplified kinetics
    ids_list: list-id list including all the ids in kinetics, reactants and products
    species_in_kinetic_law: list-species in the kinetics

    Returns
    -------
    True or False
    """
    flag_numerator = False
    flag_denominator = False
    flag = False
    eq = self._numeratorDenominator(kinetics_sim, ids_list)
    numerator = eq[0]
    denominator = eq[1]
    species = species_in_kinetic_law[0]

    if "+" not in numerator and "-" not in numerator:
      if species in numerator:
        if ("pow(" in numerator and "-1)" not in numerator) or "**" in numerator:
          flag_numerator = True
    if "+" in denominator:
      terms = denominator.split("+")
      term1 = terms[0]
      term2 = terms[1]
      if species in term1 and species not in term2:
        if ("pow(" in term1 and "-1)" not in term1) or "**" in term1:
          flag_denominator = True
      if species in term2 and species not in term1:
        if ("pow(" in term2 and "-1)" not in term2) or "**" in term2:
          flag_denominator = True

    if flag_numerator == True and flag_denominator == True:
      flag = True

    return flag

  def _numeratorDenominator(self, kinetics_sim, ids_list):
    """
    Get the numerator and denominator of a "fraction" function.
    
    Parameters
    ----    
    kinetics_sim: string-simplified kinetics
    ids_list: list-id list including all the ids in kinetics, reactants and products

    Returns
    -------
    Type - the numerator and the denominator of the fraction
    """

    strange_func = 0 
    pre_symbols = ''
    for i in range(len(ids_list)):
      pre_symbols += ids_list[i]
      pre_symbols += ' '
    pre_symbols = pre_symbols[:-1] #remove the space at the end
    pre_symbols_comma = pre_symbols.replace(" ",",")
    stmt = "%s = symbols('%s')"%(pre_symbols_comma,pre_symbols)
    try: #sometimes there is "invalid syntax error"
      exec(stmt,globals())
    except: 
      strange_func = 1
    
    try: #check if there is strange func (i.e. delay) in kinetic law; 
      #sometimes ids_list is not enough for all the ids in kinetics
      eq_stat = "kinetics_eq = " + kinetics_sim
      exec(eq_stat,globals())
    except:
      strange_func = 1

    eq = ['', '']
    if strange_func == 0:
      try: 
        numerator = str(kinetics_eq.as_numer_denom()[0])
        denominator = str(kinetics_eq.as_numer_denom()[1])
        eq[0] = numerator
        eq[1] = denominator
      except:
        pass

    return eq


  def _isPolynomial(self, kinetics_sim, ids_list):
    """
    Check if a function is polynomial.
    
    Parameters
    ----    
    kinetics_sim: string-simplified kinetics
    ids_list: list-id list including all the ids in kinetics, reactants and products

    Returns
    -------
    Type - True or False
    """

    strange_func = 0 #check if there are strang functions (i.e. delay) in kinetics
    pre_symbols = ''
    for i in range(len(ids_list)):
      pre_symbols += ids_list[i]
      pre_symbols += ' '
    pre_symbols = pre_symbols[:-1] #remove the space at the end
    pre_symbols_comma = pre_symbols.replace(" ",",")
    stmt = "%s = symbols('%s')"%(pre_symbols_comma,pre_symbols)
    try: #sometimes there is "invalid syntax error"
      exec(stmt,globals())
    except: 
      strange_func = 1
    
    try: #check if there is strange func (i.e. delay) in kinetic law
      eq_stat = "kinetics_eq = " + kinetics_sim
      exec(eq_stat,globals())
    except:
      strange_func = 1

    polynomial_flag = False
    if strange_func == 0:
      try:
        polynomial_flag = kinetics_eq.is_polynomial()
      except:
        pass
    return polynomial_flag


  @staticmethod
  def _expandFormula(expansion, function_definitions,
        num_recursion=0):
    """
    Expands the kinetics formula, replacing function definitions
    with their body.

    Parameters
    ----------
    expansion: str
        expansion of the kinetic law
    function_definitions: list-FunctionDefinition
    num_recursion: int
    
    Returns
    -------
    str
    """
    if num_recursion > MAX_RECURSION:
      return expansion
    done = True
    for fd in function_definitions:
      # Find the function calls
      calls = re.findall(r'{}\(.*?\)'.format(fd.id), expansion)
      if len(calls) == 0:
        continue
      done = False
      for call in calls:
        # Find argument call. Ex: '(a, b)'
        call_arguments = re.findall(r'\(.*?\)', call)[0]
        call_arguments = call_arguments.strip()
        call_arguments = call_arguments[1:-1]  # Eliminate parentheses
        arguments = call_arguments.split(',')
        arguments = [a.strip() for a in arguments]
        body = str(fd.body)
        for formal_arg, call_arg in zip(fd.argument_names, arguments):
          body = body.replace(formal_arg, call_arg)
        expansion = expansion.replace(call, body)
    if not done:
      return KineticLaw._expandFormula(expansion, function_definitions,
          num_recursion=num_recursion+1)
    return expansion
 

  def _getSymbols(self):
    """
    Finds the parameters and species names for the
    kinetics law. Exposing this information requires
    a recursive search of the parse tree for the
    kinetics expression.
    :return list-str:
    """

    global cur_depth
    MAX_DEPTH = 20
    cur_depth = 0
    def augment(ast_node, result):
      global cur_depth
      cur_depth += 1
      #This step is necessary to avoid memory leaking
      if cur_depth > MAX_DEPTH:
        #self.reaction.id is the reason for "try"
        raise exceptions.BadKineticsMath(self.reaction.id)
      for idx in range(ast_node.getNumChildren()):
        child_node = ast_node.getChild(idx)
        if child_node.getName() is None:
          additions = augment(child_node, result)
          result.extend(additions)
        else:
          if child_node.isFunction():
            additions = augment(child_node, result)
            result.extend(additions)
          else:
            result.append(child_node.getName())
      return result

    ast_node = self.libsbml_kinetics.getMath()
    if ast_node.getName() is None:
      result = []
    else:
      result = [ast_node.getName()]
    return augment(ast_node, result)
