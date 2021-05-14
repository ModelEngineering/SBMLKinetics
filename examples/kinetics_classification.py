
"""
This script is to do kinetic classification.
Make sure that you have setup your PYTHONPATH environment
variable as described in the github repository.
"""

# Import the required files
from sympy.core import parameters
from src.common.simple_sbml import SimpleSBML
import src.common.simple_sbml as simple_sbml
import src.common.constants as cn

import numpy as np
import collections #use set to compare two lists
import os

import sympy
from sympy import *

from libsbml import *
import tesbml # access functions in SBML
import re # Extract substrings between brackets
import time
start_time = time.time()

initial = 11

iterator = simple_sbml.modelIterator(initial=initial, final= 12)

#do statistics for different types of reactions and non-classified reactions
rxn_num = 0        #total number of reactions deals
rxn_zero_num = 0   #reaction number for zero order reaction
rxn_hill_num = 0   #reaction number for kinetics with hill terms
rxn_exp_num = 0    #reaction number for kinetics with exponential terms
rxn_no_rct_num = 0 #reaction number w/o reactants and not above
rxn_uni_num = 0    #reaction number for uni-directional mass reaction
rxn_uni_mod_num = 0#reaction number for uni-term with moderator
rxn_bi_num = 0     #reaction number for bi-directional mass reaction
rxn_bi_mod_num = 0 #reaction number for bi-terms with moderator
rxn_mm_num = 0     #reaction number for Michaelis-Menton kinetics
rxn_mm_cat_num = 0 #reaction number for Michaelis-Menton-catalyzed kinetics
rxn_mm_rev_num = 0 #reaction number for Michaelis-Menton-reversible kinetics
rxn_no_prd_num = 0 #reaction number w/o products and not above

file = open("classification.txt", "w+")
file.write("SBML id \tReaction id \tReaction \tKinetic law \tTypes of kinetics \n")

file_mol_stat = open("statistics_per_model.txt", "w+")
file_mol_stat.write("SBMLid \tReaction# \tzeroth \tHill \tExp \tno_rct \tuni \tuni_mod \tbi \
                               \tbi_mod \tMM \tMM_cat \tMM_rev \tno_prd \tNA \n")  
for idx, item in enumerate(iterator):
  if item is None:
    file_num = initial+idx
    print("File %d has an error." % (file_num))
  else:
    name = item.filename
    #print(name)
    print(name[10:])
    
    # Create an SBML model. We'll use the model
    # data/
    path = os.path.join(cn.PROJECT_DIR, "data")
    path = os.path.join(path, name)
    simple = item.model # Create a model

    model = simple.model
    # If there are functions in the sbml file, expand the functions to kinetic law first
    if (model.getNumFunctionDefinitions() != 0):
      func_id_list = [] #function id list
      var_list = []     #list of variables in the function
      func_list = []    #the equation of the function
      for f in range(model.getNumFunctionDefinitions()):
        var_list_per_func = []
        fd = model.getFunctionDefinition(f)
        func_id_list.append(fd.getId())
        for arg in range(fd.getNumArguments()):
          var_list_per_func.append(fd.getArgument(arg).getName())
        var_list.append(var_list_per_func)
        func_list.append(tesbml.formulaToL3String(fd.getBody()))

      for reaction in simple.reactions:    
        #check if there are functions in the certain reaction
        possible_func = []
        for i in range(model.getNumFunctionDefinitions()):
          if func_id_list[i] in reaction.kinetic_law.formula:
            possible_func.append(func_id_list[i])
  
        func_trial = 0
        #there are funcs in funcs, so give five trials
        while len(possible_func) > 0 and func_trial < 5: 
          #to obtain the function name, removing the error when there is only part of the 
          #function name in the kinetics    
          func_trial += 1
          func_len = len(possible_func[0])
          KL_func_id = possible_func[0]
          for i in range(1,len(possible_func)):
            if len(possible_func[i]) > func_len:
              KL_func_id = possible_func[i]
              func_len = len(possible_func[i])
          #do the function expansion for the certain kinetics
          for i in range(model.getNumFunctionDefinitions()):
            if func_id_list[i] == KL_func_id:
              #search for the "func_id(variables)" and obtain the list of variables
              bw_bracket = re.findall(r'{}\(.*?\)'.format(KL_func_id), reaction.kinetic_law.formula)
              if len(bw_bracket) > 0:
                bw_bracket = re.findall(r'\(.*?\)', bw_bracket[0])
                if len(bw_bracket) > 0:
                  KL_var_list_per_func = bw_bracket[0][1:-1].split(',')
                  #exchange any variables shown differently from the functions in the 
                  #front of the sbml file.
                  if (len(var_list[i]) == len(KL_var_list_per_func)):
                    temp = func_list[i].replace(var_list[i][0],KL_var_list_per_func[0])
                    for j in range(1,len(var_list[i])):
                      temp = temp.replace(var_list[i][j],KL_var_list_per_func[j][1:])
                    KL_func_full = KL_func_id+bw_bracket[0]
                    reaction.kinetic_law.formula = reaction.kinetic_law.formula.replace(KL_func_full,temp)
          possible_func = []
          for i in range(model.getNumFunctionDefinitions()):
            if func_id_list[i] in reaction.kinetic_law.formula:
              possible_func.append(func_id_list[i])        

    #do the statistics per model
    rxn_num_permol = len(simple.reactions)
    if rxn_num_permol != 0:
      file_mol_stat.write("%s \t" % name[10:])
      file_mol_stat.write("%s \t" % rxn_num_permol)
      rxn_zero_num_permol = 0   
      rxn_hill_num_permol = 0  
      rxn_exp_num_permol = 0      
      rxn_no_rct_num_permol = 0
      rxn_uni_num_permol = 0    
      rxn_uni_mod_num_permol = 0
      rxn_bi_num_permol = 0     
      rxn_bi_mod_num_permol = 0 
      rxn_mm_num_permol = 0     
      rxn_mm_cat_num_permol = 0 
      rxn_mm_rev_num_permol = 0 
      rxn_no_prd_num_permol = 0
      for reaction in simple.reactions:
        #change for the exponential from ^ to **
        reaction.kinetic_law.formula = reaction.kinetic_law.formula.replace('^','**')
        file.write("%s \t" % name[10:])
        rxn_num += 1
        file.write("%s \t" % reaction.getId())
        reactant_list = []
        product_list = []

        reactant_stg = " + ".join(
          [r.getSpecies() for r in reaction.reactants])
        reactant_list.append([r.getSpecies() for r in reaction.reactants])
        product_stg = " + ".join(
          [p.getSpecies() for p in reaction.products])
        product_list.append([p.getSpecies() for p in reaction.products])

        print("%s -> %s; %s" % (
          reactant_stg, product_stg,
          reaction.kinetic_law.formula))


        file.write("%s -> %s \t" % (
          reactant_stg, product_stg))

        file.write("%s \t" % (reaction.kinetic_law.formula))

        species_num = model.getNumSpecies()
        parameter_num = model.getNumParameters()

        species_list = []
        parameter_list = []
        for i in range(species_num):
          species = model.getSpecies(i)
          species_id = species.getId()
          species_list.append(species_id)

        for i in range(parameter_num):
          parameter = model.getParameter(i)
          parameter_id =  parameter.getId()
          parameter_list.append(parameter_id)

        ids_list = list(dict.fromkeys(reaction.kinetic_law.symbols))

        #type: zeroth order
        #classification rule: if there are no species in the kinetics
        if all(s not in ids_list for s in species_list):
          print("zeroth order")
          file.write("zeroth order")
          rxn_zero_num_permol += 1

        else:
          kinetics = reaction.kinetic_law.formula
          #tpye: kinetics with hill terms
          #classification rule: if there is pow() or ** inside the kinetics, 
          #except the pow(,-1) case as the possible Michaelis–Menten kinetics.
          if "pow(" in kinetics and "-1)" not in kinetics:
            print("kinetics with hill terms")
            file.write("kinetics with hill terms")
            rxn_hill_num_permol += 1            
          elif "**" in kinetics:        
            print("kinetics with hill terms")
            file.write("kinetics with hill terms")
            rxn_hill_num_permol += 1
          #type: kinetics with exponential terms
          #classification rule: if there is exp() inside the kinetics  
          elif "exp(" in kinetics:
            print("kinetics with exponential terms")
            file.write("kinetics with exponential terms")
            rxn_exp_num_permol += 1 
          #type: no reactants
          #classifcation rule: if there are no reactants
          elif len(reactant_list[0]) == 0:
            print("no reactants and not above")
            file.write("no reactants and not above")
            rxn_no_rct_num_permol += 1
          else:
            strange_func = 0 #check if there are strang functions (i.e. delay) in kinetics
            species_in_kinetic_law = []
            parameters_in_kinetic_law = []
            others_in_kinetic_law = []

            for i in range(len(ids_list)):
              if ids_list[i] in species_list:
                species_in_kinetic_law.append(ids_list[i])
              elif ids_list[i] in parameter_list:
                parameters_in_kinetic_law.append(ids_list[i])
              else:
                others_in_kinetic_law.append(ids_list[i])
            
            parameters_in_kinetic_law = parameters_in_kinetic_law + others_in_kinetic_law
            # print("species")
            # print(species_in_kinetic_law)
            # print("parameters")
            # print(parameters_in_kinetic_law)

            kinetics = reaction.kinetic_law.formula  
            #type: uni-term including uni-directional mass reaction and uni-term with moderator
            #classification rule: there is only * inside the kinetics without /,+,-.
            #for uni-directional mass reaction: the species inside the kinetics are only reactants     
            if "/" not in kinetics and "+" not in kinetics and "-" not in kinetics:
              if collections.Counter(species_in_kinetic_law) == collections.Counter(reactant_list[0]):
                print("uni-directional mass reaction")
                file.write("uni-directional mass reaction")
                rxn_uni_num_permol += 1
              else:
                print("uni-term with moderator")
                file.write("uni-term with moderator")
                rxn_uni_mod_num_permol += 1 
            #type: bi-term including bi-directional mass reaction and bi-term with moderator
            #classification rule: there is only *,- inside the kinetics without /,+.
            #for the bi-directional mass reaction: the first term before - includes all the reactants,
            #while the second term after - includes all the products.
            elif "/" not in kinetics and "+" not in kinetics and "-" in kinetics:
              terms = kinetics.split("-")
              if len(terms) == 2:
                term1 = terms[0]
                term2 = terms[1]
                if collections.Counter(species_in_kinetic_law) == collections.Counter(reactant_list[0]+product_list[0]):
                  if all(ele in term1 for ele in reactant_list[0]) and all(ele in term2 for ele in product_list[0]):
                    print("bi-directional mass reaction")
                    file.write("bi-directional mass reaction")
                    rxn_bi_num_permol += 1
                  else:
                    print("bi-terms with moderator")
                    file.write("bi-terms with moderator")
                    rxn_bi_mod_num_permol += 1 
                else:
                  print("bi-terms with moderator")
                  file.write("bi-terms with moderator")
                  rxn_bi_mod_num_permol += 1
          
            else:
              if len(reactant_list[0]) != 0:
                ids_list += reactant_list[0] # some rcts/prds also needs symbols definition

              if len(product_list[0]) != 0:
                ids_list += product_list[0]
            
              ids_list = list(dict.fromkeys(ids_list))

              pre_symbols = ''
              for i in range(len(ids_list)):
                pre_symbols += ids_list[i]
                pre_symbols += ' '
              pre_symbols = pre_symbols[:-1] #remove the space at the end
              pre_symbols_comma = pre_symbols.replace(" ",",")
              stmt = "%s = symbols('%s')"%(pre_symbols_comma,pre_symbols)
              try: #sometimes there is "invalid syntax error"
                exec(stmt)
              except: 
                strange_func = 1

              try: #check if there is strange func (i.e. delay) in kinetic law
                kinetics = reaction.kinetic_law.formula
                expr_stat = "expr = " + kinetics
                exec(expr_stat)
              except:
                strange_func = 1
              if strange_func == 0:
                #double check hill/uni/bi after simplifying and expanding the kinetics
                #following the rules stated above
                kinetics = reaction.kinetic_law.formula
                try:
                  kinetics = str(simplify(kinetics))
                except:
                  kinetics = kinetics
                
                try: #for hill terms
                  kinetics = str(expand(kinetics))
                except:
                  kinetics = kinetics
 
                if "**" in kinetics:        
                  print("kinetics with hill terms")
                  file.write("kinetics with hill terms")
                  rxn_hill_num_permol += 1
                elif "/" not in kinetics and "+" not in kinetics and "-" not in kinetics:
                  if collections.Counter(species_in_kinetic_law) == collections.Counter(reactant_list[0]):
                    print("uni-directional mass reaction")
                    file.write("uni-directional mass reaction")
                    rxn_uni_num_permol += 1
                  else:
                    print("uni-term with moderator")
                    file.write("uni-term with moderator")
                    rxn_uni_mod_num_permol += 1 

                elif "/" not in kinetics and "+" not in kinetics and "-" in kinetics:
                  terms = kinetics.split("-")
                  if len(terms) == 2:
                    term1 = terms[0]
                    term2 = terms[1]
                    if collections.Counter(species_in_kinetic_law) == collections.Counter(reactant_list[0]+product_list[0]):
                      if all(ele in term1 for ele in reactant_list[0]) and all(ele in term2 for ele in product_list[0]):
                        print("bi-directional mass reaction")
                        file.write("bi-directional mass reaction")
                        rxn_bi_num_permol += 1
                      else:
                        print("bi-terms with moderator")
                        file.write("bi-terms with moderator")
                        rxn_bi_mod_num_permol += 1 
                    else:
                      print("bi-terms with moderator")
                      file.write("bi-terms with moderator")
                      rxn_bi_mod_num_permol += 1
                else:
                  #Michaelis–Menten kinetics(inreversible)
                  #classification rule:assuming there are one/two/three parameters in the numerator,
                  #use "simplify" equals to
                  flag_mm = 0
                  if (len(reactant_list[0])==1 and len(species_in_kinetic_law)==1):
                    if (len(parameters_in_kinetic_law) >= 2):                    
                      for j in range(len(parameters_in_kinetic_law)):
                        for k in range(len(parameters_in_kinetic_law)):
                          for l in range(len(parameters_in_kinetic_law)):
                            for m in range(len(parameters_in_kinetic_law)):
                              # assuming there is one parameter in the numerator
                              if k != j:
                                pre_n = reactant_list[0][0]
                                pre_d = ' ( '
                                pre_d += reactant_list[0][0] 
                                pre_n += ' * '
                                pre_d += ' + '
                                pre_n += parameters_in_kinetic_law[j]
                                pre_d += parameters_in_kinetic_law[k]
                                pre_d += ' ) '
                                pre = pre_n
                                pre += ' / '
                                pre += pre_d           
                                expr1_stat = "expr1 =" + pre
                                exec(expr1_stat)
                                if simplify(expr1) == simplify(expr):
                                  flag_mm = 1
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
                                    exec(expr1_stat)
                                    if simplify(expr1) == simplify(expr):
                                      flag_mm = 1
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
                                        exec(expr1_stat)
                                        if simplify(expr1) == simplify(expr):
                                          flag_mm = 1
                                          break
                            if flag_mm == 1:
                              break
                          if flag_mm == 1:
                            break
                        if flag_mm == 1:
                          break
                  if flag_mm == 1:                     
                    print("Michaelis-Menten Kinetics")
                    file.write("Michaelis-Menten Kinetics")
                    rxn_mm_num_permol += 1  
                  
                  else:
                    #Michaelis–Menten kinetics(catalyzed)
                    #classification rule:assuming there are no/one/two parameters in the numerator,
                    #use "simplify" equals to
                    flag_mm_cat = 0
                    if (len(reactant_list[0])==1 and len(species_in_kinetic_law)==2):
                      if (len(parameters_in_kinetic_law) != 0):                    
                        for j in range(len(parameters_in_kinetic_law)):
                          for k in range(len(parameters_in_kinetic_law)):
                            for l in range(len(parameters_in_kinetic_law)):
                              #no parameter in the numerator
                              pre_n = reactant_list[0][0]
                              cat = [item for item in species_in_kinetic_law if item not in reactant_list[0]][0]
                              pre_n += ' * '
                              pre_n += cat
                              pre_d = ' ( '
                              pre_d += reactant_list[0][0] 
                              pre_d += ' + '
                              pre_d += parameters_in_kinetic_law[k]
                              pre_d += ' ) '
                              pre = pre_n
                              pre += ' / '
                              pre += pre_d           
                              expr1_stat = "expr1 =" + pre
                              exec(expr1_stat)
                              if simplify(expr1) == simplify(expr):
                                flag_mm_cat = 1
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
                                  exec(expr1_stat)
                                  if simplify(expr1) == simplify(expr):
                                    flag_mm_cat = 1
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
                                      exec(expr1_stat)
                                      if simplify(expr1) == simplify(expr):
                                        flag_mm_cat = 1
                                        break
                            if flag_mm_cat == 1:
                              break
                          if flag_mm_cat == 1:
                            break        
                    if flag_mm_cat == 1:                     
                      print("Michaelis-Menten Kinetics-catalyzed")
                      file.write("Michaelis-Menten Kinetics-catalyzed")
                      rxn_mm_cat_num_permol += 1 
                    else:
                      #Michaelis–Menten kinetics(reversible)
                      #classification rule:assuming there no/one/two catalyzation terms,
                      #only considering one parameter in the numerator,
                      #use "simplify" equals to
                      flag_mm_rev = 0
                      if (len(reactant_list[0])==1 and len(product_list[0])==1):
                        if (len(species_in_kinetic_law) == 2): #w/o cat
                          if (len(parameters_in_kinetic_law) >= 2):                    
                            for j in range(len(parameters_in_kinetic_law)):
                              for k in range(len(parameters_in_kinetic_law)):
                                for l in range(len(parameters_in_kinetic_law)):
                                  for m in range(len(parameters_in_kinetic_law)):
                                  # assuming there is one parameter in the numerator
                                    if k != j:
                                      rct_n = reactant_list[0][0]
                                      rct_d = ' ( '
                                      rct_d += reactant_list[0][0] 
                                      rct_n += ' * '
                                      rct_d += ' + '
                                      rct_n += parameters_in_kinetic_law[j]
                                      rct_d += parameters_in_kinetic_law[k]
                                      rct_d += ' ) '
                                      rct = rct_n
                                      rct += ' / '
                                      rct += rct_d 
                                      if l != m:
                                        prd_n = product_list[0][0]
                                        prd_d = ' ( '
                                        prd_d += product_list[0][0] 
                                        prd_n += ' * '
                                        prd_d += ' + '
                                        prd_n += parameters_in_kinetic_law[l]
                                        prd_d += parameters_in_kinetic_law[m]
                                        prd_d += ' ) '
                                        prd = prd_n
                                        prd += ' / '
                                        prd += prd_d
                                        expr1_stat = "expr1 =" + rct + '-' + prd
                                        exec(expr1_stat)
                                        if simplify(expr1) == simplify(expr):
                                          flag_mm_rev = 1 
                                          break
                                  if flag_mm_rev == 1:
                                    break
                                if flag_mm_rev == 1:
                                  break
                              if flag_mm_rev == 1:
                                break      
                        else: # with cat 
                          reactant_product_list = reactant_list[0] + product_list[0]                
                          if (len(parameters_in_kinetic_law) >= 2):                    
                            for j in range(len(parameters_in_kinetic_law)):
                              for k in range(len(parameters_in_kinetic_law)):
                                for l in range(len(parameters_in_kinetic_law)):
                                  for m in range(len(parameters_in_kinetic_law)):
                                  # assuming there is one parameter in the numerator
                                    if k != j:
                                      rct_n = reactant_list[0][0]
                                      rct_d = ' ( '
                                      rct_d += reactant_list[0][0] 
                                      rct_n += ' * '
                                      rct_d += ' + '
                                      rct_n += parameters_in_kinetic_law[j]
                                      rct_d += parameters_in_kinetic_law[k]
                                      rct_d += ' ) '
                                      rct = rct_n
                                      rct += ' / '
                                      rct += rct_d 
                                      if l != m:
                                        prd_n = product_list[0][0]
                                        prd_d = ' ( '
                                        prd_d += product_list[0][0] 
                                        prd_n += ' * '
                                        prd_d += ' + '
                                        prd_n += parameters_in_kinetic_law[l]
                                        prd_d += parameters_in_kinetic_law[m]
                                        prd_d += ' ) '
                                        prd = prd_n
                                        prd += ' / '
                                        prd += prd_d
                                        cat = [item for item in species_in_kinetic_law if item not in reactant_product_list]
                                        if len(cat) == 1: # one term w cat 
                                          expr1_stat = "expr1 =" + cat[0] + '*' + rct  + '-' + prd
                                          exec(expr1_stat)
                                          if simplify(expr1) == simplify(expr):
                                            flag_mm_rev = 1
                                            break
                                          expr1_stat = "expr1 =" + rct + '-' + cat[0] + '*' + prd
                                          exec(expr1_stat)
                                          if simplify(expr1) == simplify(expr):
                                            flag_mm_rev = 1
                                            break
                                        elif len(cat) == 2:# two terms with cat
                                          expr1_stat = "expr1 =" + cat[0] + '*' + rct + '-' + cat[1] + '*' + prd
                                          exec(expr1_stat)
                                          if simplify(expr1) == simplify(expr):
                                            flag_mm_rev = 1
                                            break
                                          expr1_stat = "expr1 =" + cat[1] + '*' + rct + '-' + cat[0] + '*' + prd
                                          exec(expr1_stat)
                                          if simplify(expr1) == simplify(expr):
                                            flag_mm_rev = 1
                                            break
                                  if flag_mm_rev == 1:
                                    break
                                if flag_mm_rev == 1:
                                  break
                              if flag_mm_rev == 1:
                                break
                      if flag_mm_rev == 1:                     
                        print("Michaelis-Menten Kinetics-reversible")
                        file.write("Michaelis-Menten Kinetics-reversible")
                        rxn_mm_rev_num_permol += 1

                      elif len(product_list[0]) == 0:
                        #type: no products and not above
                        #classification rule: if there are no products
                        print("no products and not above")
                        file.write("no products and not above")
                        rxn_no_prd_num_permol += 1

                      # else:
                      #  print("\n\n")
                      #  print("not classified")
                      #  print("%s -> %s; %s" % (
                      #      reactant_stg, product_stg,
                      #      reaction.kinetic_law.formula))
                      #  print("\n\n")
              else:
                #make up the no products case which have been ignored due to strange functions.
                if len(product_list[0]) == 0:
                  print("no products and not above")
                  file.write("no products and not above")
                  rxn_no_prd_num_permol += 1

        file.write("\n")
      rxn_zero_num += rxn_zero_num_permol
      rxn_hill_num += rxn_hill_num_permol
      rxn_exp_num += rxn_exp_num_permol
      rxn_no_rct_num += rxn_no_rct_num_permol
      rxn_uni_num += rxn_uni_num_permol
      rxn_uni_mod_num += rxn_uni_mod_num_permol
      rxn_bi_num += rxn_bi_num_permol
      rxn_bi_mod_num += rxn_bi_mod_num_permol 
      rxn_mm_num += rxn_mm_num_permol
      rxn_mm_cat_num += rxn_mm_cat_num_permol
      rxn_mm_rev_num += rxn_mm_rev_num_permol
      rxn_no_prd_num += rxn_no_prd_num_permol
      rxn_non_num_permol = rxn_num_permol - rxn_zero_num_permol - rxn_hill_num_permol \
                        - rxn_exp_num_permol - rxn_no_rct_num_permol \
                        - rxn_uni_num_permol - rxn_uni_mod_num_permol - rxn_bi_num_permol \
                        - rxn_bi_mod_num_permol - rxn_mm_num_permol - rxn_mm_cat_num_permol \
                        - rxn_mm_rev_num_permol - rxn_no_prd_num_permol
      file_mol_stat.write("%f \t" % float(rxn_zero_num_permol/rxn_num_permol))
      file_mol_stat.write("%f \t" % float(rxn_hill_num_permol/rxn_num_permol))
      file_mol_stat.write("%f \t" % float(rxn_exp_num_permol/rxn_num_permol))
      file_mol_stat.write("%f \t" % float(rxn_no_rct_num_permol/rxn_num_permol))
      file_mol_stat.write("%f \t" % float(rxn_uni_num_permol/rxn_num_permol))
      file_mol_stat.write("%f \t" % float(rxn_uni_mod_num_permol/rxn_num_permol))
      file_mol_stat.write("%f \t" % float(rxn_bi_num_permol/rxn_num_permol))
      file_mol_stat.write("%f \t" % float(rxn_bi_mod_num_permol/rxn_num_permol))
      file_mol_stat.write("%f \t" % float(rxn_mm_num_permol/rxn_num_permol))
      file_mol_stat.write("%f \t" % float(rxn_mm_cat_num_permol/rxn_num_permol))
      file_mol_stat.write("%f \t" % float(rxn_mm_rev_num_permol/rxn_num_permol))
      file_mol_stat.write("%f \t" % float(rxn_no_prd_num_permol/rxn_num_permol))
      file_mol_stat.write("%f \n" % float(rxn_non_num_permol/rxn_num_permol))

file_mol_stat.close()
file.close()

if(rxn_num != 0):
  print("\n\n")
  print("brief classified reaction statistics:")
  file_gen_stat = open("general_statistics.txt", "w+")
  file_gen_stat.write("brief classified reaction statistics:\n")
  rxn_non_num = rxn_num - rxn_zero_num - rxn_hill_num - rxn_exp_num \
                        - rxn_no_rct_num - rxn_uni_num - rxn_uni_mod_num - rxn_bi_num - rxn_bi_mod_num \
                        - rxn_mm_num - rxn_mm_cat_num - rxn_mm_rev_num \
                        - rxn_no_prd_num

  file_gen_stat.write("reaction number: %d \n" % rxn_num)
  file_gen_stat.write("zeroth order: %f \n" % float(rxn_zero_num/rxn_num))
  file_gen_stat.write("Kinetics with Hill terms: %f \n" % float(rxn_hill_num/rxn_num))
  file_gen_stat.write("Kinetics with exponential terms: %f \n" % float(rxn_exp_num/rxn_num))
  file_gen_stat.write("no reactants and not above: %f \n" % float(rxn_no_rct_num/rxn_num))
  file_gen_stat.write("uni-directional mass reaction: %f \n" % float(rxn_uni_num/rxn_num))
  file_gen_stat.write("uni-term with moderator: %f \n" % float(rxn_uni_mod_num/rxn_num))
  file_gen_stat.write("bi-directional mass reaction: %f \n" % float(rxn_bi_num/rxn_num))
  file_gen_stat.write("bi-terms with moderator: %f \n" % float(rxn_bi_mod_num/rxn_num))
  file_gen_stat.write("Michaelis-Menten Kinetics: %f \n" % float(rxn_mm_num/rxn_num))
  file_gen_stat.write("Michaelis-Menten Kinetics-catalyzed: %f \n" % float(rxn_mm_cat_num/rxn_num))
  file_gen_stat.write("Michaelis-Menten Kinetics-reversible: %f \n" % float(rxn_mm_rev_num/rxn_num))
  file_gen_stat.write("no products and not above: %f \n" % float(rxn_no_prd_num/rxn_num))
  file_gen_stat.write("not classified reactions: %f \n" % float(rxn_non_num/rxn_num))
  file_gen_stat.close()

  # print("reaction number:", rxn_num)
  # print("zeroth order:", float(rxn_zero_num/rxn_num))
  # print("Kinetics with Hill terms:", float(rxn_hill_num/rxn_num))
  # print("Kinetics with exponential terms:", float(rxn_exp_num/rxn_num))
  # print("no reactants and not above:", float(rxn_no_rct_num/rxn_num))
  # print("uni-directional mass reaction:", float(rxn_uni_num/rxn_num))
  # print("uni-term with moderator:", float(rxn_uni_mod_num/rxn_num))
  # print("bi-directional mass reaction:", float(rxn_bi_num/rxn_num))
  # print("bi-terms with moderator:", float(rxn_bi_mod_num/rxn_num))
  # print("Michaelis-Menten Kinetics:", float(rxn_mm_num/rxn_num))
  # print("Michaelis-Menten Kinetics-catalyzed:", float(rxn_mm_cat_num/rxn_num))
  # print("Michaelis-Menten Kinetics-reversible:", float(rxn_mm_rev_num/rxn_num))
  # print("no products and not above:", float(rxn_no_prd_num/rxn_num))
  print("not classified reactions:", float(rxn_non_num/rxn_num))

else:
  print("There are no reactions.")


print("--- %s seconds ---" % (time.time() - start_time))
