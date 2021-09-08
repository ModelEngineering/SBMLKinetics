
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

import pandas as pd

# Column names
CLASSIFICATION_NAME = "classification name"


def main(initial_model_indx, final_model_indx): 
  """
  Process the classification of kinetics for BioModel dataset
  
  input
  -------
  initial_model_indx: int-the intial BioModel to process
  final_model_indx: int-the final BioModel to process
  
  Returns
  -------
  df_classification: DataFrame-the kinetics classification for each reaction
  df_gen_stat: DataFrame-the general statistics for all the BioModels
  df_mol_stat: DataFrame-the statistics for each BioModel
  """
  iterator = simple_sbml.modelIterator(initial=initial_model_indx, final=final_model_indx)

  #do statistics for different types of reactions and non-classified reactions
  rxn_num = 0        #total number of reactions deals
  #all the lists are following the same order of kinetics classifications
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

  df_classification = pd.DataFrame(columns = ['SBML id', 'Reaction id',  'Classifications', 'Reaction', \
            'Kinetic law', 'Zeroth order', 'Kinetics with Hill terms', \
            'No products', 'No reactants', 'Single reactant', 'Multiple reactants', \
            'Uni-directional mass reaction', 'Uni-term with moderator', \
            'Bi-directional mass reaction', 'Bi-terms with moderator', \
            'Michaelis-Menten kinetics', 'Michaelis-Menten kinetics-catalyzed', 'NA'])
  df_classification_row = 0
  df_mol_stat = pd.DataFrame(columns=['SBMLid', 'Reaction#', 'Zeroth', 'Hill', 'no_prd', 'No_rct',\
   'Sig_rct', 'Mul_rct', 'uni', 'uni_mod', 'bi', 'bi_mod', 'MM', 'MM_cat', 'NA'])
  df_mol_stat_row = 0

  for idx, item in enumerate(iterator):
    if item is None:
      file_num = initial_model_indx +idx
      print("File %d has an error." % (file_num))
    else:
      name = item.filename

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
        list_mol_stat = []
        list_mol_stat.append(name)
        list_mol_stat.append(rxn_num_permol)
        rxn_classification_num_permol = [0]*(num_type_classification+1)

        for reaction in simple.reactions:          
          flag_classification = [0]*num_type_classification
          flag_non = 1
          classification_list = []
          reaction.kinetic_law.mkSymbolExpression(simple.function_definitions)
          list_file = []
          list_file.append(name)
          list_file.append(reaction.getId())
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

          classification_cp = [
            reaction.kinetic_law.isZerothOrder(**kwargs),
            reaction.kinetic_law.isHillTerms(**kwargs),
            reaction.kinetic_law.isNoPrds(**kwargs),
            reaction.kinetic_law.isNoRcts(**kwargs),
            reaction.kinetic_law.isSingleRct(**kwargs),
            reaction.kinetic_law.isMulRcts(**kwargs),
            reaction.kinetic_law.isUNDR(**kwargs),
            reaction.kinetic_law.isUNMO(**kwargs),
            reaction.kinetic_law.isBIDR(**kwargs),
            reaction.kinetic_law.isBIMO(**kwargs),
            reaction.kinetic_law.isMM(**kwargs),
            reaction.kinetic_law.isMMcat(**kwargs),
          ]

          for i in range(num_type_classification):
            if classification_cp[i]:
              rxn_classification_num_permol[i] += 1
              flag_classification[i] = 1
              flag_non = 0
              classification_list.append(types_simplified_name[i])

          classification_str = ','.join([str(e) for e in classification_list])
          list_file.append(classification_str)
          list_file.append(str(reaction_str))
          list_file.append(reaction.kinetic_law.expanded_formula)

          for i in range(num_type_classification): 
            if flag_classification[i] == 1:
              list_file.append('x')
            else:
              list_file.append('')

          if flag_non == 1:
            rxn_classification_num_permol[num_type_classification] += 1
            list_file.append('x')
          else:
            list_file.append('')

          df_classification.loc[df_classification_row] = list_file
          df_classification_row += 1

        for i in range(num_type_classification+1):
          rxn_classification_num[i] += rxn_classification_num_permol[i] 

        rxn_num += rxn_num_permol

        for i in range(num_type_classification):
          list_mol_stat.append(float(rxn_classification_num_permol[i]/rxn_num_permol))
        list_mol_stat.append(float(rxn_classification_num_permol[num_type_classification]/rxn_num_permol))
        df_mol_stat.loc[df_mol_stat_row] = list_mol_stat
        df_mol_stat_row += 1

  # This part is the same as the printed part in main section
  if(rxn_num != 0):
    df_gen_stat = pd.DataFrame(columns=[CLASSIFICATION_NAME, 'Percentage'])
    for i in range(num_type_classification+1):
      df_gen_stat.loc[i] = [types_name[i], float(rxn_classification_num[i]/rxn_num)]

  return (df_classification, df_gen_stat, df_mol_stat)


if __name__ == '__main__':
  start_time = time.time()
  initial_model_indx = 0
  final_model_indx = 2
  (df_classification, df_gen_stat, df_mol_stat) = main(initial_model_indx, final_model_indx)
  rxn_num = len(df_classification)
  
  SBML_id_list = []
  for i in range(len(df_classification)):
    SBML_id = df_classification.iloc[i]['SBML id']
    if SBML_id not in SBML_id_list:
      print(SBML_id)
      SBML_id_list.append(SBML_id)
    print(df_classification.iloc[i]['Reaction'] + ";" + df_classification.iloc[i]['Kinetic law'])
    print(df_classification.iloc[i]['Classifications'])

  if(rxn_num != 0):
    print("\n\n")
    print("brief classified reaction statistics:")
    print("Reaction number:", rxn_num)
    for i in range(len(df_gen_stat)):
      print(df_gen_stat.iloc[i]['Classification Names'] + ":" + str(df_gen_stat.iloc[i]['Percentage']))
  else:
    print("There are no reactions.")

  df_classification.to_csv("classification.csv", index=False)
  df_gen_stat.to_csv("general_statistics.csv", index=False)
  df_mol_stat.to_csv("statistics_per_model.csv", index=False)
  print("--- %s seconds ---" % (time.time() - start_time))
