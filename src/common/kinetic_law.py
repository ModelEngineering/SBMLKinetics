'''Provides Information on SBML Kinetics Laws'''


from src.common import constants as cn
from src.common import util
from src.common import exceptions
from src.common import msgs
from sympy import *

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
    self.formula = self.libsbml_kinetics.getFormula()
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
    Creates a string that can be processed by pySym.
    
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

  def isMulRcts(self, **kwargs):
    """
    Tests whether the reaction belongs to the type of Multiple Reactants
    Multiple reactants classification rule: if there are multiple reactants
    
    Parameters
    -------
    **kwargs: dictionary-keyword arguments
    reactant_list: list-reactants of the reaction
    
    Returns
    -------
    True or False
    """
    reactant_list = kwargs["reactant_list"] 
    return self._numOfRcts(reactant_list) > 1


  def isUNDR(self, **kwargs):
    """
    Tests whether the reaction belongs to the type of uni-directional mass reaction
    Uni-directional mass reaction classification rule: 
    1) There is only * inside the kinetics without /,+,-.
    2) The species inside the kinetics are only reactants
    
    Parameters
    -------
    **kwargs: dictionary-keyword arguments
    reactant_list: list-reactants of the reaction
    kinetics: string-kinetics
    kinetics_sim: string-simplified kinetics
    species_in_kinetic_law: list-species in the kinetics
    
    Returns
    -------
    True or False
    """
    reactant_list = kwargs["reactant_list"]
    kinetics = kwargs["kinetics"]
    kinetics_sim = kwargs["kinetics_sim"]
    species_in_kinetic_law = kwargs["species_in_kinetic_law"]

    flag = True
    # if self._numSpeciesInKinetics(species_in_kinetic_law) == 0:
    #   flag = False
    if self._isSingleProductOfTerms(kinetics, kinetics_sim) == False:
      flag = False
    if self._SpecsInKineticsAllRcts(species_in_kinetic_law, reactant_list) == False:
      flag = False

    return flag

  def isUNMO(self, **kwargs):
    """
    Tests whether the reaction belongs to the type of uni-term with moderator
    Uni-term with moderator classification rule: 
    1) There is only * inside the kinetics without /,+,-.
    2) The species inside the kinetics are not only reactants
    
    Parameters
    -------
    **kwargs: dictionary-keyword arguments
    reactant_list: list-reactants of the reaction
    kinetics: string-kinetics
    kinetics_sim: string-simplified kinetics
    species_in_kinetic_law: list-species in the kinetics
    
    Returns
    -------
    True or False
    """
    reactant_list = kwargs["reactant_list"]
    kinetics = kwargs["kinetics"]
    kinetics_sim = kwargs["kinetics_sim"]
    species_in_kinetic_law = kwargs["species_in_kinetic_law"]


    flag = True
    if self._numSpeciesInKinetics(species_in_kinetic_law) == 0:
      flag = False
    if self._isSingleProductOfTerms(kinetics, kinetics_sim) == False:
      flag = False
    if self._SpecsInKineticsAllRcts(species_in_kinetic_law, reactant_list) == True:
      flag = False
    
    return flag

  def isBIDR(self, **kwargs):
    """
    Tests whether the reaction belongs to the type of bi-directional mass reaction
    Bi-directional mass reactionclassification rule: 
    1) There is only *,- inside the kinetics without /,+.
    2) The first term before - includes all the reactants,
       while the second term after - includes all the products. 
      (Is there a better and more accurate way for this?)
    
    Parameters
    -------
    **kwargs: dictionary-keyword arguments
    reactant_list: list-reactants of the reaction
    product_list: list-products of the reaction
    kinetics: string-kinetics
    kinetics_sim: string-simplified kinetics
    species_in_kinetic_law: list-species in the kinetics
    
    Returns
    -------
    True or False
    """

    reactant_list = kwargs["reactant_list"]
    product_list = kwargs["product_list"]
    kinetics = kwargs["kinetics"]
    kinetics_sim = kwargs["kinetics_sim"]
    species_in_kinetic_law = kwargs["species_in_kinetic_law"]

    flag = True
    # if self._numSpeciesInKinetics(species_in_kinetic_law) == 0:
    #   flag = False
    if self._isDiffOfTwoProductsOfTerms(kinetics, kinetics_sim) == False:
      flag = False
    if self._ProductOfTermsWithAllRctsOrPrds(kinetics, kinetics_sim, species_in_kinetic_law, reactant_list, product_list) == False:
      flag = False
    
    return flag

  def isBIMO(self, **kwargs):
    """
    Tests whether the reaction belongs to the type of bi-terms with moderator
    Bi-terms with moderator classification rule: 
    1) There is only *,- inside the kinetics without /,+.
    2) The first term before - does not include all the reactants,
       while the second term after - does not include all the products. 
      (Is there a better and more accurate way for this?)
    
    Parameters
    -------
    **kwargs: dictionary-keyword arguments
    reactant_list: list-reactants of the reaction
    product_list: list-products of the reaction
    kinetics: string-kinetics
    kinetics_sim: string-simplified kinetics
    species_in_kinetic_law: list-species in the kinetics
    
    Returns
    -------
    True or False
    """
    reactant_list = kwargs["reactant_list"]
    product_list = kwargs["product_list"]
    kinetics = kwargs["kinetics"]
    kinetics_sim = kwargs["kinetics_sim"]
    species_in_kinetic_law = kwargs["species_in_kinetic_law"]

    flag = True
    if self._numSpeciesInKinetics(species_in_kinetic_law) == 0:
      flag = False
    if self._isDiffOfTwoProductsOfTerms(kinetics, kinetics_sim) == False:
      flag = False
    if self._ProductOfTermsWithAllRctsOrPrds(kinetics, kinetics_sim, species_in_kinetic_law, reactant_list, product_list) == True:
      flag = False

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
    ids_list: list-id list including all the ids in kinetics, reactants and products
    species_in_kinetic_law: list-species in the kinetics
    parameters_in_kinetic_law: list-parameters in the kinetics  
    reactant_list: list-reactants of the reaction
    
    Returns
    -------
    True or False
    """
  
    kinetics = kwargs["kinetics"]
    ids_list = kwargs["ids_list"]
    species_in_kinetic_law = kwargs["species_in_kinetic_law"]
    parameters_in_kinetic_law = kwargs["parameters_in_kinetic_law"]
    reactant_list = kwargs["reactant_list"]

    flag = False
    if self._numSpeciesInKinetics(species_in_kinetic_law) == 1 and self._numOfRcts(reactant_list) == 1:
      if self._MMSingleSpecInNumerator(kinetics, ids_list, parameters_in_kinetic_law, reactant_list) == True:
        flag = True

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
    ids_list = kwargs["ids_list"]
    species_in_kinetic_law = kwargs["species_in_kinetic_law"]
    parameters_in_kinetic_law = kwargs["parameters_in_kinetic_law"]
    reactant_list = kwargs["reactant_list"]

    flag = False
    if self._numSpeciesInKinetics(species_in_kinetic_law) == 2 and self._numOfRcts(reactant_list) == 1:
      if self._MMTwoSpecInNumerator(kinetics, ids_list, parameters_in_kinetic_law, species_in_kinetic_law, reactant_list) == True:  
        flag = True
      
    return flag

  def isHyperbolic(self, **kwargs):
    """
    Tests whether the reaction belongs to the type of Hyperbolic function
    Hyperbolic function classification rule:
    Check whether it is in the form of fraction and has species in both numerator and denominator.
    
    Parameters
    ----  
    **kwargs: dictionary-keyword arguments  
    #kinetics: string-kinetics
    kinetics_sim: string-simplified kinetics
    reactant_list: list-reactants of the reaction
    ids_list: list-id list including all the ids in kinetics, reactants and products

    Returns
    -------
    True or False
    """
    
    #kinetics = kwargs["kinetics"]
    kinetics_sim = kwargs["kinetics_sim"]
    reactant_list = kwargs["reactant_list"]
    ids_list = kwargs["ids_list"]

    eq = self._isFraction(kinetics_sim, ids_list)
    flag = False
    if len(reactant_list) > 0:
      for i in range(len(reactant_list)):
        if reactant_list[i] in eq[0] and reactant_list[i] in eq[1]:
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
    if "/" not in kinetics and "+" not in kinetics and "-" not in kinetics:
      return True
    elif "/" not in kinetics_sim and "+" not in kinetics_sim and "-" not in kinetics_sim:
      return True
    else:
      return False 

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
    if "/" not in kinetics and "+" not in kinetics and "exp(-" not in kinetics and "-" in kinetics:
      return True
    elif "/" not in kinetics_sim and "+" not in kinetics_sim and "exp(-" not in kinetics_sim and "-" in kinetics_sim:     
      return True
    else:
      return False

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
                  if simplify(expr1) == simplify(expr):
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
                      if simplify(expr1) == simplify(expr):
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
                          if simplify(expr1) == simplify(expr):
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
              if simplify(expr1) == simplify(expr):
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
                  if simplify(expr1) == simplify(expr):
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
                      if simplify(expr1) == simplify(expr):
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

  def _isFraction(self, kinetics_sim, ids_list):
    """
    Test whether the simplified kinetics is a fraction
    
    Parameters
    ----    
    kinetics_sim: string-simplified kinetics
    ids_list: list-id list including all the ids in kinetics, reactants and products

    Returns
    -------
    Type - the numerator and the denominator of the fraction
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
