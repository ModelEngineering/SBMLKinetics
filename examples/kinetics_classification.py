
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
import os

import sympy
from sympy import symbols, simplify 

from libsbml import *
import tesbml # access functions in SBML
import re # Extract substrings between brackets
import time
start_time = time.time()

initial = 1

iterator = simple_sbml.modelIterator(initial=initial, final= 2)

#do statistics for different types of reactions and non-classified reactions
rxn_num = 0        #total number of reactions deals
rxn_zero_num = 0   #reaction number for zero order reaction
rxn_no_rct_num = 0     #reaction number w/o reactants (excluding zero order reaction)
rxn_uni_num = 0    #reaction number for uni-directional mass reaction
rxn_bi_num = 0     #reaction number for bi-directional mass reaction
rxn_mm_num = 0     #reaction number for Michaelis-Menton kinetics
rxn_mm_cat_num = 0 #reaction number for Michaelis-Menton-catalyzed kinetics
rxn_mm_rev_num = 0 #reaction number for Michaelis-Menton-catalyzed kinetics

file = open("classification.txt", "w+")
file.write("SBML id \tReaction id \tReaction \tKinetic law \tTypes of kinetics \n")
  
for idx, item in enumerate(iterator):
  if item is None:
    file_num = initial+idx
    print("File %d has an error." % (file_num))
  else:
    name = item.filename
    print(name)
    #print(name[10:])
    
    # Create an SBML model. We'll use the model
    # data/
    path = os.path.join(cn.PROJECT_DIR, "data")
    path = os.path.join(path, name)
    simple = item.model # Create a model

    model = simple.model

    # try to access the info of functions,
    # If there are functions, expand the functions to kinetic law first
    if (model.getNumFunctionDefinitions() != 0):
      func_id_list = []
      var_list = []
      func_list = []
      for f in range(model.getNumFunctionDefinitions()):
        var_list_per_func = []
        fd = model.getFunctionDefinition(f)
        func_id_list.append(fd.getId())
        for arg in range(fd.getNumArguments()):
          var_list_per_func.append(fd.getArgument(arg).getName())
        var_list.append(var_list_per_func)
        func_list.append(tesbml.formulaToL3String(fd.getBody()))

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

      for reaction in simple.reactions:

        reactant_stg = " + ".join(
          [r.getSpecies() for r in reaction.reactants])
        product_stg = " + ".join(
          [p.getSpecies() for p in reaction.products])

        possible_func = []
        for i in range(model.getNumFunctionDefinitions()):
          if func_id_list[i] in reaction.kinetic_law.formula:
            possible_func.append(func_id_list[i])
        
        if len(possible_func) > 0:
          func_len = len(possible_func[0])
          KL_func_id = possible_func[0]
          for i in range(1,len(possible_func)):
            if len(possible_func[i]) > func_len:
              KL_func_id = possible_func[i]

          for i in range(model.getNumFunctionDefinitions()):
            if func_id_list[i] == KL_func_id:
              bw_bracket = re.findall(r'\(.*?\)', reaction.kinetic_law.formula)
              if (len(bw_bracket) > 0):
                KL_var_list_per_func = bw_bracket[0][1:-1].split(',')
                if (len(var_list[i]) == len(KL_var_list_per_func)):
                  temp = func_list[i].replace(var_list[i][0],KL_var_list_per_func[0])
                  for j in range(1,len(var_list[i])):
                    temp = temp.replace(var_list[i][j],KL_var_list_per_func[j][1:])
                  KL_func_full = KL_func_id+bw_bracket[0]
                  reaction.kinetic_law.formula = reaction.kinetic_law.formula.replace(KL_func_full,temp)

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

      # print("%s -> %s; %s" % (
      #   reactant_stg, product_stg,
      #   reaction.kinetic_law.formula))


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

      if all([s not in ids_list for s in species_list]):
        print("zeroth order")
        file.write("zeroth order")
        rxn_zero_num += 1

      else:
        if len(reactant_list[0]) == 0:
          print("no reactants and not zero order")
          file.write("no reactants and not zero order")
          rxn_no_rct_num += 1
        else:
          strange_func = 0 #check if there is strang func (i.e. delay) in kinetic law
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

          reactant_times = ""
          product_times = ""
          if len(reactant_list[0]) != 0:
            ids_list += reactant_list[0] # some rcts/prds also needs definition
            reactant_times += reactant_list[0][0]
            for i in range(1,len(reactant_list[0])):
              reactant_times += " * "
              reactant_times += reactant_list[0][i]

          if len(product_list[0]) != 0:
            ids_list += product_list[0]
            product_times += product_list[0][0]
            for i in range(1, len(product_list[0])):
              product_times += " * "
              product_times += product_list[0][i]
        
          ids_list = list(dict.fromkeys(ids_list))

          pre_symbols = ''
          for i in range(len(ids_list)):
            pre_symbols += ids_list[i]
            pre_symbols += ' '
          pre_symbols = pre_symbols[:-1] #remove the space at the end

          pre_symbols_comma = pre_symbols.replace(" ",",")
          stmt = "%s = symbols('%s')"%(pre_symbols_comma,pre_symbols)
          try:
            exec(stmt)
          except: 
            strange_func = 1

          expr_stat = reaction.kinetic_law.formula
          expr_stat = "expr = " + reaction.kinetic_law.formula
          try: #check if there is strang func (i.e. delay) in kinetic law
            exec(expr_stat)
          except:
            strange_func = 1

          if strange_func == 0:
            flag_uni = 0
            if(len(reactant_times) != 0):
              if (len(parameters_in_kinetic_law) != 0):
                for j in range(len(parameters_in_kinetic_law)):
                  pre_rct = reactant_times
                  pre_rct += ' * '
                  pre_rct += parameters_in_kinetic_law[j]
                  expr1_stat = "expr1 =" + pre_rct
                  exec(expr1_stat)
                  if simplify(expr1) == simplify(expr):
                    flag_uni = 1
                  for k in range(len(parameters_in_kinetic_law)):
                    if len(parameters_in_kinetic_law) >= 2 :
                      if j!=k :
                        pre_rct += ' * '
                        pre_rct += parameters_in_kinetic_law[k]
                        expr1_stat = "expr1 =" + pre_rct
                        exec(expr1_stat)
                        if simplify(expr1) == simplify(expr):
                          flag_uni = 1
                      for l in range(len(parameters_in_kinetic_law)):
                        if len(parameters_in_kinetic_law) >= 3:
                          if l != j and l != k:
                            pre_rct += ' * '
                            pre_rct += parameters_in_kinetic_law[l]
                            expr1_stat = "expr1 =" + pre_rct
                            exec(expr1_stat)
                            if simplify(expr1) == simplify(expr):
                              flag_uni = 1
                              break
                    if flag_uni == 1:
                      break
                  if flag_uni == 1:
                    break
                  
                      
            if flag_uni == 1:
              print("uni-directional mass reaction")
              file.write("uni-directional mass reaction")
              rxn_uni_num += 1

            else:
              flag_bi = 0


              if(len(reactant_times) != 0 and len(product_times) != 0):
                if (len(parameters_in_kinetic_law) != 0):
                  for j in range(len(parameters_in_kinetic_law)):
                    for k in range(len(parameters_in_kinetic_law)):
                      for l in range(len(parameters_in_kinetic_law)):
                        for m in range(len(parameters_in_kinetic_law)):
                          for n in range(len(parameters_in_kinetic_law)):
                            pre_rct = reactant_times
                            pre_rct += ' * '
                            pre_rct += parameters_in_kinetic_law[j]
                          
                            pre_pdt = product_times
                            pre_pdt += ' * '
                            pre_pdt += parameters_in_kinetic_law[k] 

                            pre_law = pre_rct 
                            pre_law += ' - '
                            pre_law += pre_pdt
                          
                            expr1_stat = "expr1 =" + pre_law
                            exec(expr1_stat)
                            if simplify(expr1) == simplify(expr):
                              flag_bi = 1
                              break 
                            if len(parameters_in_kinetic_law) >= 2:
                              if l != j and m != k:
                                pre_rct += ' * '
                                pre_rct += parameters_in_kinetic_law[l]

                                pre_pdt += ' * '
                                pre_pdt += parameters_in_kinetic_law[m] 
                                pre_law = pre_rct 
                                pre_law += ' - '
                                pre_law += pre_pdt
                            
                                expr1_stat = "expr1 =" + pre_law
                                exec(expr1_stat)
                                if simplify(expr1) == simplify(expr):
                                  flag_bi = 1
                                  break
                                if n!=l and n!=j:  
                                  pre_rct += ' * '
                                  pre_rct += parameters_in_kinetic_law[n] 
                                  pre_law = pre_rct 
                                  pre_law += ' - '
                                  pre_law += pre_pdt
                                  
                                  expr1_stat = "expr1 =" + pre_law
                                  exec(expr1_stat)
                                  if simplify(expr1) == simplify(expr):
                                    flag_bi = 1
                                    break
                                if n!=m and n!=k:
                                  pre_rct = reactant_times
                                  pre_rct += ' * '
                                  pre_rct += parameters_in_kinetic_law[j]
                                  pre_rct += ' * '
                                  pre_rct += parameters_in_kinetic_law[l]
                                  pre_pdt = product_times
                                  pre_pdt += ' * '
                                  pre_pdt += parameters_in_kinetic_law[k]
                                  pre_pdt += ' * ( '
                                  pre_pdt += parameters_in_kinetic_law[m]
                                  pre_pdt += ' + '
                                  pre_pdt += parameters_in_kinetic_law[n]
                                  pre_pdt += ' ) ' 
                                  pre_law = pre_rct 
                                  pre_law += ' - '
                                  pre_law += pre_pdt
                                  
                                  expr1_stat = "expr1 =" + pre_law
                                  exec(expr1_stat)
                                  if simplify(expr1) == simplify(expr):
                                    flag_bi = 1
                                    break
                          if flag_bi == 1:
                            break
                        if flag_bi == 1:
                          break
                      if flag_bi == 1:
                        break
                    if flag_bi == 1:
                      break                  
              if flag_bi == 1:
                print("bi-directional mass reaction")
                file.write("bi-directional mass reaction")
                rxn_bi_num += 1
          
              else:
                flag_mm = 0
                if (len(reactant_list[0])==1):
                  if (len(parameters_in_kinetic_law) >= 2):                    
                    for j in range(len(parameters_in_kinetic_law)):
                      for k in range(len(parameters_in_kinetic_law)):
                        for l in range(len(parameters_in_kinetic_law)):
                          for m in range(len(parameters_in_kinetic_law)):
                            # assuming there are one parameter in the numerator
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
                                  expr1_stat = "expr1 =" + pr
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
                  rxn_mm_num += 1  
                
                else:
                  flag_mm_cat = 0
                  if (len(reactant_list[0])==1 and len(species_in_kinetic_law)==2):
                    if (len(parameters_in_kinetic_law) != 0):                    
                      for j in range(len(parameters_in_kinetic_law)):
                        for k in range(len(parameters_in_kinetic_law)):
                          for l in range(len(parameters_in_kinetic_law)):
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
                            # assuming there are one parameter in the numerator
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
                    rxn_mm_cat_num += 1 
                  else:
                    flag_mm_rev = 0
                    if (len(reactant_list[0])==1 and len(product_list[0])==1):
                      if (len(species_in_kinetic_law) == 2):
                        if (len(parameters_in_kinetic_law) >= 2):                    
                          for j in range(len(parameters_in_kinetic_law)):
                            for k in range(len(parameters_in_kinetic_law)):
                              for l in range(len(parameters_in_kinetic_law)):
                                for m in range(len(parameters_in_kinetic_law)):
                                # assuming there are one parameter in the numerator
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
                      else:  
                        reactant_product_list = reactant_list[0] + product_list[0]                
                        if (len(parameters_in_kinetic_law) >= 2):                    
                          for j in range(len(parameters_in_kinetic_law)):
                            for k in range(len(parameters_in_kinetic_law)):
                              for l in range(len(parameters_in_kinetic_law)):
                                for m in range(len(parameters_in_kinetic_law)):
                                # assuming there are one parameter in the numerator
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
                                      if len(cat) == 1:
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
                                      elif len(cat) == 2:
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
                      rxn_mm_rev_num += 1

                    #else:
                    #  print("\n\n")
                    #  print("not classified")
                    #  print("%s -> %s; %s" % (
                    #      reactant_stg, product_stg,
                    #      reaction.kinetic_law.formula))
                    #  print("\n\n")

      file.write("\n")      

if(rxn_num != 0):
  print("\n\n")
  print("brief classified reaction statistics:")
  rxn_non_num = rxn_num - rxn_zero_num - rxn_no_rct_num - rxn_uni_num - rxn_bi_num -rxn_mm_num-rxn_mm_cat_num-rxn_mm_rev_num
  print("reaction number:", rxn_num)
  print("zeroth order:", float(rxn_zero_num/rxn_num))
  print("no reactants and not zero order:", float(rxn_no_rct_num/rxn_num))
  print("uni-directional mass reaction:", float(rxn_uni_num/rxn_num))
  print("bi-directional mass reaction:", float(rxn_bi_num/rxn_num))
  print("Michaelis-Menten Kinetics:", float(rxn_mm_num/rxn_num))
  print("Michaelis-Menten Kinetics-catalyzed:", float(rxn_mm_cat_num/rxn_num))
  print("Michaelis-Menten Kinetics-reversible:", float(rxn_mm_rev_num/rxn_num))
  print("not classified reactions:", float(rxn_non_num/rxn_num))

else:
  print("There are no reactions.")

file.close()

print("--- %s seconds ---" % (time.time() - start_time))
