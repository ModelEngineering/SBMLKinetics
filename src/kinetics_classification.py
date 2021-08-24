
"""
This script is to do kinetic classification.
Make sure that you have setup your PYTHONPATH environment
variable as described in the github repository.
"""

# Import the required files
import types
from sympy.core import parameters
from src.common.simple_sbml import SimpleSBML
import src.common.simple_sbml as simple_sbml
import src.common.constants as cn

import numpy as np
import os

from sympy import *
from libsbml import * # access functions in SBML
import time



def main(initial_model_indx, final_model_indx): 
  iterator = simple_sbml.modelIterator(initial=initial_model_indx, final=final_model_indx)

  #do statistics for different types of reactions and non-classified reactions
  rxn_num = 0        #total number of reactions deals
  types_name = ["Zeroth Order", "Kinetics with Hill terms", "No products", \
    "No reactants", "Single reactant", "Multiple reactants", \
    "Uni-directional mass reaction", "Uni-term with moderator", \
    "Bi-directional mass reaction", "Bi-terms with moderator", \
    "Michaelis-Menten Kinetics", "Michaelis-Menten Kinetics-catalyzed", \
    "Not classified reactions"]
  types_simplified_name = ["ZERO", "HILL", "P=0", "R=0", "R=1", "R>1", \
    "UNDR", "UNMO", "BIDR", "BIMO", "MM", "MMCAT"]
  num_type_classification = len(types_name) - 1  #types_name includes not classified cases
  rxn_classification_num = [0]*(num_type_classification+1)

  file = open("classification.txt", "w+")
  file.write("SBML id \tReaction id  \tClassifications \tReaction \tKinetic law \
            \tZeroth order \tKinetics with Hill terms \
            \tNo products \tNo reactants \tSingle reactant \tMultiple reactants \
            \tUni-directional mass reaction \tUni-term with moderator \
            \tBi-directional mass reaction \tBi-terms with moderator \
            \tMichaelis-Menten kinetics \tMichaelis-Menten kinetics-catalyzed \
            \tNA\n")
  file_mol_stat = open("statistics_per_model.txt", "w+")
  file_mol_stat.write("SBMLid \tReaction# \tZeroth \tHill \tno_prd \tNo_rct \tSig_rct \tMul_rct \
                        \tuni \tuni_mod \tbi \tbi_mod \tMM \tMM_cat \tNA \n") 

  for idx, item in enumerate(iterator):
    if item is None:
      file_num = initial_model_indx +idx
      print("File %d has an error." % (file_num))
    else:
      name = item.filename
      print(name)
      #print(name[10:])

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
        rxn_classification_num_permol = [0]*(num_type_classification+1)

        for reaction in simple.reactions:          
          flag_classification = [0]*num_type_classification
          flag_non = 1
          classification_list = []
          reaction.kinetic_law.mkSymbolExpression(simple.function_definitions)
          file.write("%s \t" % name)
          file.write("%s \t" % reaction.getId())
          reactant_list = [r.getSpecies() for r in reaction.reactants]
          product_list = [p.getSpecies() for p in reaction.products]

          #print("reactant_list")
          #print(reactant_list)
          #print("product_list")
          #print(product_list)

          reactant_stg = " + ".join(
            [r.getSpecies() for r in reaction.reactants])
          product_stg = " + ".join(
            [p.getSpecies() for p in reaction.products])

          print(str(reaction))
          reaction_str = reactant_stg + "->" + product_stg

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

          #only for MM and MMcat
          if len(reactant_list) != 0:
            ids_list += reactant_list # some rcts/prds also needs symbols definition
          if len(product_list) != 0:
            ids_list += product_list
          ids_list = list(dict.fromkeys(ids_list))

          # Define the keyword arguments
          kwargs = {"kinetics": kinetics, "kinetics_sim": kinetics_sim, \
            "reactant_list": reactant_list, "product_list": product_list, \
            "species_in_kinetic_law": species_in_kinetic_law, "parameters_in_kinetic_law": parameters_in_kinetic_law, \
            "ids_list": ids_list}

          ZERO_cp = reaction.kinetic_law.isZerothOrder(**kwargs)
          HILL_cp = reaction.kinetic_law.isHillTerms(**kwargs)
          NO_P_cp = reaction.kinetic_law.isNoPrds(**kwargs)
          NO_R_cp = reaction.kinetic_law.isNoRcts(**kwargs)
          SIG_R_cp = reaction.kinetic_law.isSingleRct(**kwargs)
          MUL_R_cp = reaction.kinetic_law.isMulRcts(**kwargs)
          UNDR_cp = reaction.kinetic_law.isUNDR(**kwargs)
          UNMO_cp = reaction.kinetic_law.isUNMO(**kwargs)
          BIDR_cp = reaction.kinetic_law.isBIDR(**kwargs)
          BIMO_cp = reaction.kinetic_law.isBIMO(**kwargs)
          MM_cp = reaction.kinetic_law.isMM(**kwargs)
          MMCAT_cp = reaction.kinetic_law.isMMcat(**kwargs)
          classification_func = [ZERO_cp, HILL_cp, NO_P_cp, NO_R_cp, SIG_R_cp, MUL_R_cp, \
                                UNDR_cp, UNMO_cp, BIDR_cp, BIMO_cp, MM_cp, MMCAT_cp]

          for i in range(num_type_classification):
            if classification_func[i]:
              rxn_classification_num_permol[i] += 1
              flag_classification[i] = 1
              flag_non = 0
              classification_list.append(types_simplified_name[i])

          classification_str = ','.join([str(elem) for elem in classification_list])
          file.write(classification_str)
          file.write("\t")

          file.write(str(reaction_str))
          file.write("\t")
          file.write("%s \t" % (reaction.kinetic_law.expanded_formula))

          for i in range(num_type_classification): 
            if flag_classification[i] == 1:
              print(types_name[i])
              file.write("x \t")
            else:
              file.write("\t")

          if flag_non == 1:
            rxn_classification_num_permol[num_type_classification] += 1
            file.write("x \n")
          else:
            file.write("\n") 

        for i in range(num_type_classification+1):
          rxn_classification_num[i] += rxn_classification_num_permol[i] 

        rxn_num += rxn_num_permol

        for i in range(num_type_classification):
          file_mol_stat.write("%f \t" % float(rxn_classification_num_permol[i]/rxn_num_permol))
        file_mol_stat.write("%f \n" % float(rxn_classification_num_permol[num_type_classification]/rxn_num_permol))

  file_mol_stat.close()
  file.close()
  if(rxn_num != 0):
    file_gen_stat = open("general_statistics.txt", "w+")
    file_gen_stat.write("brief classified reaction statistics:\n")
    file_gen_stat.write("Reaction number: %d \n" % rxn_num)
    for i in range(num_type_classification+1):
      file_gen_stat.write(types_name[i] + ": %f \n" % float(rxn_classification_num[i]/rxn_num))
    file_gen_stat.close()

  return (types_name, rxn_classification_num, rxn_num)


if __name__ == '__main__':
  start_time = time.time()
  initial_model_indx = 0
  final_model_indx = 10
  (types_name, rxn_classification_num, rxn_num) = main(initial_model_indx, final_model_indx)
  if(rxn_num != 0):
    print("\n\n")
    print("brief classified reaction statistics:")
    print("Reaction number:", rxn_num)
    for i in range(len(types_name)):
      print(types_name[i] + ":", float(rxn_classification_num[i]/rxn_num))
  else:
    print("There are no reactions.")
  print("--- %s seconds ---" % (time.time() - start_time))
