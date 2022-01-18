
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
import math

# Column names
SBMLID = "SBMLid"
REACTIONID = "Reaction id"
CLASSIFICATIONS = 'Classifications'
REACTION = 'Reaction'
KINETICLAW = 'kinetic law'
ZEROTH = 'Zeroth order'
POWER = 'Kinetics with power terms'
NOPRD = 'No products'
NORCT = 'No reactants'
SIGRCT = 'Single reactant'
MULRCT = 'Multiple reactants'
UNI = 'Uni-directional mass reaction'
UNIMOD = 'Uni-term with moderator'
BI = 'Bi-directional mass reaction'
BIMOD = 'Bi-terms with moderator'
MM = 'Michaelis-Menten kinetics'
MMCAT = 'Michaelis-Menten kinetics-catalyzed'
FR = "Fraction"
PL = "Polynomial"
NA = 'NA'
PERCENTAGE = 'Percentage'
PERCENTAGE_SDER = 'Percentage standard error'
PERCENTAGE_PER_MODEL = 'Percentage per model'
PERCENTAGE_PER_MODEL_SDER = 'Percentage per model standard error'
RXN_NUM = 'Reaction number'
BIOMOL_NUM = 'Biomodel number'
# COLUMN_NAME_df_classification = [SBMLID, REACTIONID, CLASSIFICATIONS, REACTION, KINETICLAW,
#                 ZEROTH, POWER, NOPRD, NORCT, SIGRCT, MULRCT, UNI, UNIMOD, BI, BIMOD, MM, MMCAT, FR, PL, NA]

COLUMN_NAME_df_classification = [SBMLID, REACTIONID, CLASSIFICATIONS, REACTION, KINETICLAW,
                ZEROTH, UNI, UNIMOD, BI, BIMOD, FR, NA]

# COLUMN_NAME_df_gen_stat = [CLASSIFICATIONS, PERCENTAGE, PERCENTAGE_SDER, \
#   PERCENTAGE_PER_MODEL, PERCENTAGE_PER_MODEL_SDER, RXN_NUM, BIOMOL_NUM]
COLUMN_NAME_df_gen_stat = [CLASSIFICATIONS, PERCENTAGE, \
  PERCENTAGE_PER_MODEL, PERCENTAGE_PER_MODEL_SDER, RXN_NUM, BIOMOL_NUM]
# COLUMN_NAME_df_mol_stat = [SBMLID, RXN_NUM, ZEROTH, POWER, NOPRD, NORCT,\
#   SIGRCT, MULRCT, UNI, UNIMOD, BI, BIMOD, MM, MMCAT, FR, PL, NA]

COLUMN_NAME_df_mol_stat = [SBMLID, RXN_NUM, ZEROTH, UNI, UNIMOD, BI, BIMOD, FR, NA]


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
  df_gen_stat_PR: DataFrame-the general statistics for all the BioModels per number of rcts and prds
  biomodel_non_count: Int-number of biomodels with non-classified rxns
  """
  iterator = simple_sbml.modelIterator(initial=initial_model_indx, final=final_model_indx)

  #do statistics for different types of reactions and non-classified reactions
  rxn_num = 0        #total number of reactions deals
  rxn_num_PR = [0]*16 #total number of reactions with certain prds and rcts (4prds*4rcts)
  #all the lists are following the same order of kinetics classifications
  # types_name = ["Zeroth Order", "Kinetics with power terms", "No products", \
  #   "No reactants", "Single reactant", "Multiple reactants", \
  #   "Uni-directional mass reaction", "Uni-term with moderator", \
  #   "Bi-directional mass reaction", "Bi-terms with moderator", \
  #   "Michaelis-Menten Kinetics", "Michaelis-Menten Kinetics-catalyzed", \
  #   "Fraction function", "Polynomial function"\
  #   "Not classified reactions"]
  types_name = ["Zeroth Order", \
    "Uni-directional mass reaction", "Uni-term with moderator", \
    "Bi-directional mass reaction", "Bi-terms with moderator", \
    "Fraction function", \
    "Not classified reactions"]
  # types_simplified_name = ["ZERO", "POWER", "P=0", "R=0", "R=1", "R>1", \
  #   "UNDR", "UNMO", "BIDR", "BIMO", "MM", "MMCAT", "FR", "PL"]
  types_simplified_name = ["ZERO", \
    "UNDR", "UNMO", "BIDR", "BIMO", "FR"]
  num_type_classification = len(types_name) -1  #types_name includes not classified cases
  rxn_classification_num = [0]*(num_type_classification+1) #total number of classified cases for each type
  rxn_classification_num_PR = np.zeros((16, (num_type_classification+1))) #16 rows for 4prds*4rcts

  df_classification = pd.DataFrame(columns = COLUMN_NAME_df_classification)
  df_mol_stat = pd.DataFrame(columns = COLUMN_NAME_df_mol_stat)
  df_mol_stat_PR = {} #set of dataframes to save mol stat for each PR
  for i in range(16):
    df_mol_stat_PR[i] = pd.DataFrame(columns = COLUMN_NAME_df_mol_stat)
  
  biomodel_non_count = 0
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
      rxn_num_permol_PR = [0]*16 #rxn numbers per model for each PR
      if rxn_num_permol != 0:
        flag_biomodel_non = 0
        mol_stat_row_dct = {k:[] for k in COLUMN_NAME_df_mol_stat}
        mol_stat_row_dct[SBMLID].append(name)
        mol_stat_row_dct[RXN_NUM].append(rxn_num_permol)

        rxn_classification_num_permol = [0]*(num_type_classification+1)
        rxn_classification_num_permol_PR = np.zeros((16,(num_type_classification+1)))

        for reaction in simple.reactions:          
          flag_classification = [0]*num_type_classification
          flag_non = 1
          classification_list = []
          reaction.kinetic_law.mkSymbolExpression(simple.function_definitions)   
          classification_row_dct = {k:[] for k in COLUMN_NAME_df_classification}
          classification_row_dct[SBMLID].append(name)
          classification_row_dct[REACTIONID].append(reaction.getId())
          reactant_list = [r.getSpecies() for r in reaction.reactants]
          product_list = [p.getSpecies() for p in reaction.products]

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
         
          #print("kinetics_sim:", kinetics_sim)

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
          # print("species_in_kinetic_law:", species_in_kinetic_law)
          # print("parameters")
          # print(parameters_in_kinetic_law)

          #only for MM, MMcat and FR
          if len(reactant_list) != 0:
            ids_list += reactant_list # some rcts/prds also needs symbols definition
          if len(product_list) != 0:
            ids_list += product_list
          ids_list = list(dict.fromkeys(ids_list))

          # print("reactant_list:", reactant_list)
          # print("ids_list:", ids_list)

          #Define the keyword arguments
          kwargs = {"kinetics": kinetics, "kinetics_sim": kinetics_sim, \
            "reactant_list": reactant_list, "product_list": product_list, \
            "species_in_kinetic_law": species_in_kinetic_law, "parameters_in_kinetic_law": parameters_in_kinetic_law, \
            "ids_list": ids_list}

          classification_cp = [
            reaction.kinetic_law.isZerothOrder(**kwargs),
            # reaction.kinetic_law.isPowerTerms(**kwargs),
            reaction.kinetic_law.isUNDR(**kwargs),
            reaction.kinetic_law.isUNMO(**kwargs),
            reaction.kinetic_law.isBIDR(**kwargs),
            reaction.kinetic_law.isBIMO(**kwargs),
            #reaction.kinetic_law.isMM(**kwargs),
            #reaction.kinetic_law.isMMcat(**kwargs),
            reaction.kinetic_law.isFraction(**kwargs),
            #reaction.kinetic_law.isPolynomial(**kwargs),
          ]

          classification_prds_cp = [
            reaction.kinetic_law.isNoPrds(**kwargs),
            reaction.kinetic_law.isSinglePrd(**kwargs),
            reaction.kinetic_law.isDoublePrds(**kwargs),
            reaction.kinetic_law.isMulPrds(**kwargs),
          ]

          classification_rcts_cp = [
            reaction.kinetic_law.isNoRcts(**kwargs),
            reaction.kinetic_law.isSingleRct(**kwargs),
            reaction.kinetic_law.isDoubleRcts(**kwargs),
            reaction.kinetic_law.isMulRcts(**kwargs),
          ]

          for i in range(num_type_classification):
            if classification_cp[i]:
              rxn_classification_num_permol[i] += 1
              flag_classification[i] = 1
              flag_non = 0
              classification_list.append(types_simplified_name[i])
              break #stop the loop once flag_non == 0, this applies to exclusive classification

          classification_str = ','.join([str(e) for e in classification_list])
          classification_row_dct[CLASSIFICATIONS].append(classification_str)
          classification_row_dct[REACTION].append(reaction_str)
          classification_row_dct[KINETICLAW].append(reaction.kinetic_law.expanded_formula)

          for i in range(num_type_classification): 
            if flag_classification[i] == 1:
              classification_row_dct[COLUMN_NAME_df_classification[5+i]].append('x')
            else:
              classification_row_dct[COLUMN_NAME_df_classification[5+i]].append('')

          if flag_non == 1:
            rxn_classification_num_permol[num_type_classification] += 1
            classification_row_dct[NA].append('x')
            flag_biomodel_non = 1
          else:
            classification_row_dct[NA].append('')
          
          for i in range(len(COLUMN_NAME_df_classification)):
            classification_row_dct[COLUMN_NAME_df_classification[i]] = classification_row_dct[COLUMN_NAME_df_classification[i]][0]
          df_classification = df_classification.append(classification_row_dct, ignore_index=True)


          for x in range(4): #4prds
            for y in range(4): #4rcts
              flag_non_PR = 1
              xy = x*4+y
    
              if classification_prds_cp[x] and classification_rcts_cp[y]:
                rxn_num_permol_PR[xy] += 1
                for i in range(num_type_classification):
                  if classification_cp[i]:
                    rxn_classification_num_permol_PR[xy,i] += 1
                    flag_non_PR = 0
                    break #stop the loop once, this applies to exclusive classification
                if flag_non_PR == 1:
                  rxn_classification_num_permol_PR[xy,num_type_classification] += 1 

        #for each reaction above here, for per model below here:
        for i in range(num_type_classification+1):
          rxn_classification_num[i] += rxn_classification_num_permol[i] 
        rxn_num += rxn_num_permol

        for i in range(num_type_classification):
          mol_stat_row_dct[COLUMN_NAME_df_mol_stat[2+i]].append(float(rxn_classification_num_permol[i]/rxn_num_permol))
        mol_stat_row_dct[NA].append(float(rxn_classification_num_permol[num_type_classification]/rxn_num_permol))
        
        for i in range(len(COLUMN_NAME_df_mol_stat)):
          mol_stat_row_dct[COLUMN_NAME_df_mol_stat[i]] = mol_stat_row_dct[COLUMN_NAME_df_mol_stat[i]][0]
        df_mol_stat = df_mol_stat.append(mol_stat_row_dct, ignore_index=True)
      
        if flag_biomodel_non == 1:
          biomodel_non_count += 1

        #PR:
        for xy in range(16):
          if rxn_num_permol_PR[xy]!= 0:
            mol_stat_PR_row_dct = {k:[] for k in COLUMN_NAME_df_mol_stat}
            mol_stat_PR_row_dct[SBMLID].append(name)
            
            mol_stat_PR_row_dct[RXN_NUM].append(rxn_num_permol_PR[xy])
            for i in range(num_type_classification+1):
              rxn_classification_num_PR[xy,i] += rxn_classification_num_permol_PR[xy,i] 
            rxn_num_PR[xy] += rxn_num_permol_PR[xy]

            for i in range(num_type_classification):
              mol_stat_PR_row_dct[COLUMN_NAME_df_mol_stat[2+i]].append(float(rxn_classification_num_permol_PR[xy,i]/rxn_num_permol_PR[xy]))
            mol_stat_PR_row_dct[NA].append(float(rxn_classification_num_permol_PR[xy,num_type_classification]/rxn_num_permol_PR[xy]))
            
            for i in range(len(COLUMN_NAME_df_mol_stat)):
              mol_stat_PR_row_dct[COLUMN_NAME_df_mol_stat[i]] = mol_stat_PR_row_dct[COLUMN_NAME_df_mol_stat[i]][0]
            df_mol_stat_PR[xy] = df_mol_stat_PR[xy].append(mol_stat_PR_row_dct, ignore_index=True) 
  #for all the biomodels below here
  # This part is the same as the printed part in main section
  if(rxn_num != 0):
    df_gen_stat = pd.DataFrame(columns = COLUMN_NAME_df_gen_stat)
    for i in range(num_type_classification+1):
      gen_stat_row_dct = {k:[] for k in COLUMN_NAME_df_gen_stat[0:-2]}
      gen_stat_row_dct[CLASSIFICATIONS].append(types_name[i])
      gen_stat_row_dct[PERCENTAGE].append(float(rxn_classification_num[i]/rxn_num))
      
      # do a statistics of df_mol_stat and save to df_gen_stat
      avg_value = df_mol_stat[COLUMN_NAME_df_mol_stat[i+2]].mean()
      sdv_value = df_mol_stat[COLUMN_NAME_df_mol_stat[i+2]].std()/math.sqrt(len(df_mol_stat.index))
      if math.isnan(sdv_value):
        sdv_value = 0.
      gen_stat_row_dct[PERCENTAGE_PER_MODEL].append(avg_value)
      gen_stat_row_dct[PERCENTAGE_PER_MODEL_SDER].append(sdv_value)
      for j in range(len(COLUMN_NAME_df_gen_stat)-2):
        gen_stat_row_dct[COLUMN_NAME_df_gen_stat[j]] = gen_stat_row_dct[COLUMN_NAME_df_gen_stat[j]][0]
      df_gen_stat = df_gen_stat.append(gen_stat_row_dct, ignore_index=True)
    
    df_gen_stat.at[0, RXN_NUM] = rxn_num
    df_gen_stat.at[0, BIOMOL_NUM] = len(df_mol_stat.index)


  #PR
  df_gen_stat_PR = pd.DataFrame(columns = COLUMN_NAME_df_gen_stat)
  for xy in range(16):
    row = len(df_gen_stat_PR.index)
    if(rxn_num_PR[xy] != 0):
      for i in range(num_type_classification+1):
        gen_stat_PR_row_dct = {k:[] for k in COLUMN_NAME_df_gen_stat[0:-2]}
        gen_stat_PR_row_dct[CLASSIFICATIONS].append(types_name[i])
        gen_stat_PR_row_dct[PERCENTAGE].append(float(rxn_classification_num_PR[xy,i]/rxn_num_PR[xy]))

        # do a statistics of df_mol_stat and save to df_gen_stat
        avg_value = df_mol_stat_PR[xy][COLUMN_NAME_df_mol_stat[i+2]].mean()
        sdv_value = df_mol_stat_PR[xy][COLUMN_NAME_df_mol_stat[i+2]].std()/math.sqrt(len(df_mol_stat_PR[xy].index))
        if math.isnan(sdv_value):
          sdv_value = 0
        gen_stat_PR_row_dct[PERCENTAGE_PER_MODEL].append(avg_value)
        gen_stat_PR_row_dct[PERCENTAGE_PER_MODEL_SDER].append(sdv_value)
        for j in range(len(COLUMN_NAME_df_gen_stat)-2):
          gen_stat_PR_row_dct[COLUMN_NAME_df_gen_stat[j]] = gen_stat_PR_row_dct[COLUMN_NAME_df_gen_stat[j]][0]    
        df_gen_stat_PR = df_gen_stat_PR.append(gen_stat_PR_row_dct, ignore_index=True)  
    df_gen_stat_PR.at[row, RXN_NUM] = rxn_num_PR[xy]
    df_gen_stat_PR.at[row, BIOMOL_NUM] = len(df_mol_stat_PR[xy].index)
    

  return (df_classification, df_gen_stat, df_mol_stat, df_gen_stat_PR, biomodel_non_count)


if __name__ == '__main__':
  start_time = time.time()
  initial_model_indx = 5
  final_model_indx = 6
  (df_classification, df_gen_stat, df_mol_stat, df_gen_stat_PR, biomodel_non_count) = \
    main(initial_model_indx, final_model_indx)
  rxn_num = len(df_classification)
  

  SBML_id_list = []
  for i in range(len(df_classification)):
    SBML_id = df_classification.iloc[i][SBMLID]
    if SBML_id not in SBML_id_list:
      print(SBML_id)
      SBML_id_list.append(SBML_id)
    print(df_classification.iloc[i][REACTION] + ";" + df_classification.iloc[i][KINETICLAW])
    print(df_classification.iloc[i][CLASSIFICATIONS])

  if(rxn_num != 0):
    print("\n\n")
    print("brief classified reaction statistics:")
    print("Reaction number:", rxn_num)
    for i in range(len(df_gen_stat)):
      print(df_gen_stat.iloc[i][CLASSIFICATIONS] + ":" + str(df_gen_stat.iloc[i][PERCENTAGE]))
  else:
    print("There are no reactions.")

  print("number of biomodels with some reactions not classified:", biomodel_non_count)

  # df_classification.to_csv("classification.csv", index=False)
  # df_gen_stat.to_csv("general_statistics.csv", index=False)
  # df_mol_stat.to_csv("statistics_per_model.csv", index=False)

  # Create a Pandas Excel writer using XlsxWriter as the engine.
  writer = pd.ExcelWriter('statistics_result.xlsx', engine='xlsxwriter')

  # Write each dataframe to a different worksheet.
  df_classification.to_excel(writer, sheet_name='classification')
  df_gen_stat.to_excel(writer, sheet_name='general_statistics')
  df_mol_stat.to_excel(writer, sheet_name='statistics_per_model')
  df_gen_stat_PR.to_excel(writer, sheet_name='general_statistics_PR')
  #df_mol_stat_PR.to_excel(writer, sheet_name='statistics_per_model_PR5')

  # Close the Pandas Excel writer and output the Excel file.
  writer.save()

  print("--- %s seconds ---" % (time.time() - start_time))
