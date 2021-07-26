
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
import libsbml # access functions in SBML
import re # Extract substrings between brackets
import time
start_time = time.time()

dMAX_ITERATION = 5  # Maximum number for iteration function expansions

initial = 41

iterator = simple_sbml.modelIterator(initial=initial, final=42)

#do statistics for different types of reactions and non-classified reactions
rxn_num = 0        #total number of reactions deals
rxn_zero_num = 0   #reaction number for zero order reaction
rxn_hill_num = 0   #reaction number for kinetics with hill term
rxn_no_prd_num = 0 #reaction number w/o products 
rxn_no_rct_num = 0 #reaction number w/o reactants
rxn_sig_rct_num = 0 #reaction number w single reactant
rxn_mul_rct_num = 0 #reaction number w multiple reactants
rxn_uni_num = 0    #reaction number for uni-directional mass reaction
rxn_uni_mod_num = 0#reaction number for uni-term with moderator
rxn_bi_num = 0     #reaction number for bi-directional mass reaction
rxn_bi_mod_num = 0 #reaction number for bi-terms with moderator
rxn_mm_num = 0     #reaction number for Michaelis-Menton kinetics
rxn_mm_cat_num = 0 #reaction number for Michaelis-Menton-catalyzed kinetics
rxn_non_num = 0    #reaction number does not fall in any types above

file = open("classification.txt", "w+")
file.write("SBML id \tReaction id \tReaction \tKinetic law \
           \tZeroth order \tKinetics with Hill terms \
           \tNo products \tNo reactants \tSingle reactant \tMultiple reactants \
           \tUni-directional mass reaction \tUni-term with moderator \
           \tBi-directional mass reaction \tBi-terms with moderator \
           \tMichaelis-Menten kinetics \tMichaelis-Menten kinetics-catalyzed \
           \tNA \tClassifications \n")
# file.write("SBML id \tReaction id  \tClassifications \tReaction \tKinetic law \
#            \tZeroth order \tKinetics with Hill terms \
#            \tNo products \tNo reactants \tSingle reactant \tMultiple reactants \
#            \tUni-directional mass reaction \tUni-term with moderator \
#            \tBi-directional mass reaction \tBi-terms with moderator \
#            \tMichaelis-Menten kinetics \tMichaelis-Menten kinetics-catalyzed \
#            \tNA\n")

file_mol_stat = open("statistics_per_model.txt", "w+")
file_mol_stat.write("SBMLid \tReaction# \tZeroth \tHill \tno_prd \tNo_rct \tSig_rct \tMul_rct \
                      \tuni \tuni_mod \tbi \tbi_mod \tMM \tMM_cat \tNA \n") 


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
    try:
      path = os.path.join(cn.PROJECT_DIR, "data")
    except:
      print("error")
    path = os.path.join(path, name)
    simple = item.model # Create a model

    model = simple.model
    # If there are functions in the sbml file, expand the functions to kinetic law first
    if len(simple.function_definitions) > 0:
      for reaction in simple.reactions:
        reaction.kinetic_law.expandFormula(simple.function_definitions)

    #do the statistics per model
    rxn_num_permol = len(simple.reactions)
    if rxn_num_permol != 0:
      file_mol_stat.write("%s \t" % name[10:])
      file_mol_stat.write("%s \t" % rxn_num_permol)
      rxn_zero_num_permol = 0   
      rxn_hill_num_permol = 0     
      rxn_no_prd_num_permol = 0
      rxn_no_rct_num_permol = 0
      rxn_sig_rct_num_permol = 0
      rxn_mul_rct_num_permol = 0
      rxn_uni_num_permol = 0    
      rxn_uni_mod_num_permol = 0
      rxn_bi_num_permol = 0     
      rxn_bi_mod_num_permol = 0 
      rxn_mm_num_permol = 0     
      rxn_mm_cat_num_permol = 0 
      rxn_non_num_permol = 0

      for reaction in simple.reactions:

        flag_zeroth = 0
        #flag_hill = 0
        flag_no_prd = 0
        flag_no_rct = 0
        flag_sig_rct = 0
        #flag_mul_rct = 0
        #flag_UNDR = 0
        #flag_UNMO = 0
        #flag_BIDR = 0
        #flag_BIMO = 0
        #flag_MM = 0
        #flag_MMCAT = 0
        flag_non = 1
        classification_list = []
        reaction.kinetic_law.mkSymbolExpression(simple.function_definitions)
        file.write("%s \t" % name[10:])
        #rxn_num += 1
        file.write("%s \t" % reaction.getId())
        # QUESTION: Why is this a list of lists?
        reactant_list = [[r.getSpecies() for r in reaction.reactants]]
        product_list = [[p.getSpecies() for p in reaction.products]]

        reactant_stg = " + ".join(
          [r.getSpecies() for r in reaction.reactants])
        product_stg = " + ".join(
          [p.getSpecies() for p in reaction.products])


        print(str(reaction))

        reaction_str = reactant_stg + "->" + product_stg
        file.write(str(reaction_str))
        file.write("\t")

        file.write("%s \t" % (reaction.kinetic_law.expanded_formula))

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

        kinetics = reaction.kinetic_law.expanded_formula  

        try:
          kinetics_sim = str(simplify(kinetics))
        except:
          kinetics_sim = kinetics


        ids_list = list(dict.fromkeys(reaction.kinetic_law.symbols))

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
        #print("species")
        #print(species_in_kinetic_law)
        #print("parameters")
        #print(parameters_in_kinetic_law)

    
        if reaction.kinetic_law.isZerothOrder(species_in_kinetic_law):
          print("Zeroth order")
          file.write("x \t")
          rxn_zero_num_permol += 1
          flag_non = 0
          flag_zeroth = 1
          classification_list.append("ZERO")
        else:
          file.write("\t")

        #tpye: kinetics with hill terms
        #classification rule: if there is pow() or ** inside the kinetics, 
        #except the pow(,-1) case as the possible Michaelis–Menten kinetics.
        if reaction.kinetic_law.isHillTerms(kinetics,kinetics_sim):
          print("Kinetics with Hill terms")
          file.write("x \t")
          rxn_hill_num_permol += 1
          flag_non = 0            
          classification_list.append("HILL")
        else:
          file.write("\t")


        #type: no products
        #classification rule; if there are no products
        if reaction.kinetic_law.isNoPrds(product_list):
          print("No products")
          file.write("x \t")
          rxn_no_prd_num_permol += 1
          flag_no_prd = 1
          flag_non = 0
          classification_list.append("P=0")
        else:
          file.write("\t")
          
        #type: no reactants
        #classifcation rule: if there are no reactants
        if reaction.kinetic_law.isNoRcts(reactant_list):
          print("No reactants")
          file.write("x \t")
          rxn_no_rct_num_permol += 1
          flag_no_rct = 1
          flag_non = 0
          classification_list.append("R=0")
        else:
          file.write("\t")

        if reaction.kinetic_law.isSingleRct(reactant_list):
          print("One reactant")
          file.write("x \t")
          rxn_sig_rct_num_permol += 1
          flag_sig_rct = 1
          flag_non = 0
          classification_list.append("R=1")
        else:
          file.write("\t")

        if reaction.kinetic_law.isMulRcts(reactant_list):
          print("Multiple reactants")
          file.write("x \t")
          rxn_mul_rct_num_permol += 1
          flag_non = 0
          classification_list.append("R>1")
        else:
          file.write("\t")


        #type: uni-term including uni-directional mass reaction and uni-term with moderator
        #classification rule: there is only * inside the kinetics without /,+,-.
        #for uni-directional mass reaction: the species inside the kinetics are only reactants     
        if reaction.kinetic_law.isUNDR(reactant_list,kinetics,kinetics_sim,species_in_kinetic_law):
            print("Uni-directional mass reaction")
            file.write("x \t")
            rxn_uni_num_permol += 1
            flag_non = 0
            classification_list.append("UNIDR")
        else:
            file.write("\t")

        if reaction.kinetic_law.isUNMO(reactant_list,kinetics,kinetics_sim,species_in_kinetic_law):
            print("Uni-term with moderator")
            file.write("x \t")
            rxn_uni_mod_num_permol += 1
            flag_non = 0
            classification_list.append("UNIMO")
        else: 
            file.write("\t") 

        #type: bi-term including bi-directional mass reaction and bi-term with moderator
        #classification rule: there is only *,- inside the kinetics without /,+.
        #for the bi-directional mass reaction: the first term before - includes all the reactants,
        #while the second term after - includes all the products. 
        #(Is there a better and more accurate way for this?)
        
        if reaction.kinetic_law.isBIDR(reactant_list,product_list,kinetics,kinetics_sim,species_in_kinetic_law):
          print("Bi-directional mass reaction")
          file.write("x \t")
          rxn_bi_num_permol += 1
          flag_non = 0
          classification_list.append("BIDR")
        else:
          file.write("\t")
        
        if reaction.kinetic_law.isBIMO(reactant_list,product_list,kinetics,kinetics_sim,species_in_kinetic_law):
          print("Bi-terms with moderator")
          file.write("x \t")
          rxn_bi_mod_num_permol += 1
          flag_non = 0
          classification_list.append("BIMO") 
        else:
          file.write("\t")


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
          expr_stat = "expr = " + kinetics
          exec(expr_stat)
        except:
          strange_func = 1

        if strange_func == 0:
        
          # Michaelis–Menten kinetics(inreversible)
          # classification rule:assuming there are one/two/three parameters in the numerator,
          # use "simplify" equals to
          flag_mm = 0
          if (flag_sig_rct==1 and len(species_in_kinetic_law)==1):
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
            file.write("x \t")
            rxn_mm_num_permol += 1
            flag_non = 0
            classification_list.append("MM")
          else:
            file.write("\t")
            
          # Michaelis–Menten kinetics(catalyzed)
          # classification rule:assuming there are no/one/two parameters in the numerator,
          # use "simplify" equals to
          flag_mm_cat = 0 
          if (flag_sig_rct==1 and len(species_in_kinetic_law)==2):
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
            file.write("x \t")
            rxn_mm_cat_num_permol += 1
            flag_non = 0
            classification_list.append("MMCAT")
          else:
            file.write("\t")

        else: 
          file.write("\t\t")


        # if reaction.kinetic_law.isMM(kinetics, ids_list, species_in_kinetic_law, parameters_in_kinetic_law, reactant_list):                      
        #   print("Michaelis-Menten Kinetics")
        #   file.write("x \t")
        #   rxn_mm_num_permol += 1
        #   flag_non = 0
        #   classification_list.append("MM")
        # else:
        #   file.write("\t")

        # if reaction.kinetic_law.isMMcat(kinetics, ids_list, species_in_kinetic_law, parameters_in_kinetic_law, reactant_list):                     
        #     print("Michaelis-Menten Kinetics-catalyzed")
        #     file.write("x \t")
        #     rxn_mm_cat_num_permol += 1
        #     flag_non = 0
        #     classification_list.append("MMCAT")
        # else:
        #   file.write("\t")
        
        if flag_non == 1:
          rxn_non_num_permol += 1
          file.write("x \t")
        else:
          file.write("\t") 
        
        #print(classification_list)
        classification_str = ','.join([str(elem) for elem in classification_list])
        file.write(classification_str)
        file.write("\n")

      rxn_zero_num    += rxn_zero_num_permol
      rxn_hill_num    += rxn_hill_num_permol
      rxn_no_prd_num  += rxn_no_prd_num_permol
      rxn_no_rct_num  += rxn_no_rct_num_permol
      rxn_sig_rct_num += rxn_sig_rct_num_permol 
      rxn_mul_rct_num += rxn_mul_rct_num_permol
      rxn_uni_num     += rxn_uni_num_permol
      rxn_uni_mod_num += rxn_uni_mod_num_permol
      rxn_bi_num      += rxn_bi_num_permol
      rxn_bi_mod_num  += rxn_bi_mod_num_permol
      rxn_mm_num      += rxn_mm_num_permol
      rxn_mm_cat_num  += rxn_mm_cat_num_permol
      rxn_non_num     += rxn_non_num_permol  
      rxn_num         += rxn_num_permol

      file_mol_stat.write("%f \t" % float(rxn_zero_num_permol/rxn_num_permol))
      file_mol_stat.write("%f \t" % float(rxn_hill_num_permol/rxn_num_permol))
      file_mol_stat.write("%f \t" % float(rxn_no_prd_num_permol/rxn_num_permol))
      file_mol_stat.write("%f \t" % float(rxn_no_rct_num_permol/rxn_num_permol))
      file_mol_stat.write("%f \t" % float(rxn_sig_rct_num_permol/rxn_num_permol))
      file_mol_stat.write("%f \t" % float(rxn_mul_rct_num_permol/rxn_num_permol))
      file_mol_stat.write("%f \t" % float(rxn_uni_num_permol/rxn_num_permol))
      file_mol_stat.write("%f \t" % float(rxn_uni_mod_num_permol/rxn_num_permol))
      file_mol_stat.write("%f \t" % float(rxn_bi_num_permol/rxn_num_permol))
      file_mol_stat.write("%f \t" % float(rxn_bi_mod_num_permol/rxn_num_permol))
      file_mol_stat.write("%f \t" % float(rxn_mm_num_permol/rxn_num_permol))
      file_mol_stat.write("%f \t" % float(rxn_mm_cat_num_permol/rxn_num_permol))
      file_mol_stat.write("%f \n" % float(rxn_non_num_permol/rxn_num_permol))

file_mol_stat.close()
file.close()

if(rxn_num != 0):
  print("\n\n")
  print("brief classified reaction statistics:")
  file_gen_stat = open("general_statistics.txt", "w+")
  file_gen_stat.write("brief classified reaction statistics:\n")

  file_gen_stat.write("Reaction number: %d \n" % rxn_num)
  file_gen_stat.write("Zeroth order: %f \n" % float(rxn_zero_num/rxn_num))
  file_gen_stat.write("Kinetics with Hill terms: %f \n" % float(rxn_hill_num/rxn_num))
  file_gen_stat.write("No products: %f \n" % float(rxn_no_prd_num/rxn_num))
  file_gen_stat.write("No reactants: %f \n" % float(rxn_no_rct_num/rxn_num))
  file_gen_stat.write("Single reactant: %f \n" % float(rxn_sig_rct_num/rxn_num))
  file_gen_stat.write("Multiple reactants: %f \n" % float(rxn_mul_rct_num/rxn_num))
  file_gen_stat.write("Uni-directional mass reaction: %f \n" % float(rxn_uni_num/rxn_num))
  file_gen_stat.write("Uni-term with moderator: %f \n" % float(rxn_uni_mod_num/rxn_num))
  file_gen_stat.write("Bi-directional mass reaction: %f \n" % float(rxn_bi_num/rxn_num))
  file_gen_stat.write("Bi-terms with moderator: %f \n" % float(rxn_bi_mod_num/rxn_num))
  file_gen_stat.write("Michaelis-Menten Kinetics: %f \n" % float(rxn_mm_num/rxn_num))
  file_gen_stat.write("Michaelis-Menten Kinetics-catalyzed: %f \n" % float(rxn_mm_cat_num/rxn_num))
  file_gen_stat.write("Not classified reactions: %f \n" % float(rxn_non_num/rxn_num))
  file_gen_stat.close()

  print("Reaction number:", rxn_num)
  print("Zeroth order:", float(rxn_zero_num/rxn_num))
  print("Kinetics with Hill terms:", float(rxn_hill_num/rxn_num))
  print("No products:", float(rxn_no_prd_num/rxn_num))
  print("No reactants:", float(rxn_no_rct_num/rxn_num))
  print("Single reactant:", float(rxn_sig_rct_num/rxn_num))
  print("Multiple reactants:", float(rxn_mul_rct_num/rxn_num))
  print("Uni-directional mass reaction:", float(rxn_uni_num/rxn_num))
  print("Uni-term with moderator:", float(rxn_uni_mod_num/rxn_num))
  print("Bi-directional mass reaction:", float(rxn_bi_num/rxn_num))
  print("Bi-terms with moderator:", float(rxn_bi_mod_num/rxn_num))
  print("Michaelis-Menten Kinetics:", float(rxn_mm_num/rxn_num))
  print("Michaelis-Menten Kinetics-catalyzed:", float(rxn_mm_cat_num/rxn_num))
  print("Not classified reactions:", float(rxn_non_num/rxn_num))

else:
  print("There are no reactions.")


print("--- %s seconds ---" % (time.time() - start_time))
