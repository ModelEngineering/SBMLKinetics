
# This script was written by Jin Xu and available on Github
# https://github.com/SunnyXu/SBMLKinetics
# This file includes all the functions to do the kinetics classification.


# Import the required files
import types
from zipfile import ZIP_FILECOUNT_LIMIT
from isort import file
from sympy.core import parameters
from SBMLKinetics.common.simple_sbml import SimpleSBML
import SBMLKinetics.common.simple_sbml as simple_sbml
import SBMLKinetics.common.constants as cn
import sys

import numpy as np
import os

from sympy import *
from libsbml import * # access functions in SBML
import time

import collections
import pandas as pd
import math
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib


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
UNI = 'Uni-directional mass action'
UNIMOD = 'Uni-term with moderator'
BI = 'Bi-directional mass action'
BIMOD = 'Bi-terms with moderator'
MM = 'Michaelis-Menten kinetics'
MMCAT = 'Michaelis-Menten kinetics-catalyzed'
HILL = "Hill equations"
FR = "Fraction" #kinetics in fraction format other than MM, MMCAT, HILL
PL = "Polynomial"
NA = 'NA'
PERCENTAGE = 'Percentage'
PERCENTAGE_SDER = 'Percentage standard error'
PERCENTAGE_PER_MODEL = 'Percentage per model'
PERCENTAGE_PER_MODEL_SDER = 'Percentage per model standard error'
RXN_NUM = 'Reaction number'
BIOMOL_NUM = 'Biomodel number'
COLUMN_NAME_df_classification = [SBMLID, REACTIONID, CLASSIFICATIONS, REACTION, KINETICLAW,
                ZEROTH, UNI, UNIMOD, BI, BIMOD, MM, MMCAT, HILL, FR, NA]

COLUMN_NAME_df_gen_stat = [CLASSIFICATIONS, PERCENTAGE, \
 PERCENTAGE_PER_MODEL, PERCENTAGE_PER_MODEL_SDER, RXN_NUM, BIOMOL_NUM]

COLUMN_NAME_df_mol_stat = [SBMLID, RXN_NUM, ZEROTH, UNI, UNIMOD, BI, BIMOD, MM, MMCAT, HILL, FR, NA]


def _dataSetStatistics(data_dir = cn.BIOMODELS_DIR, zip_filename = cn.BIOMODELS_ZIP_FILENAME,
initial_model_indx = 0, final_model_indx = 1000): 
  """
  Process the classification of kinetics for BioModel dataset.
  
  input
  -------
  data_dir: folder path.
  zip_filename: str-zip filfile name, e.g. "dataSetName.zip".
  initial_model_indx: int-the intial BioModel to process.
  final_model_indx: int-the final BioModel to process.
  
  Returns
  -------
  df_classification: DataFrame-the kinetics classification for each reaction.
  df_gen_stat: DataFrame-the general statistics for all the BioModels.
  df_mol_stat: DataFrame-the statistics for each BioModel.
  df_gen_stat_PR: DataFrame-the general statistics for all the BioModels per number of rcts and prds.
  biomodel_non_count: Int-number of biomodels with non-classified rxns.
  df_table_PR: DataFrame-the statistics for percentage of reactions in each type of mass transfer.
  df_table_PR_per_model-df_table_PR averagely for each model.
  """

  #do statistics for different types of reactions and non-classified reactions
  rxn_num = 0        #total number of reactions deals
  rxn_num_PR = [0]*16 #total number of reactions with certain prds and rcts (4prds*4rcts)
  #all the lists are following the same order of kinetics classifications
  types_name = ["ZERO", "UNDR", "UNMO", "BIDR", "BIMO", "MM", "MMCAT", "HILL", "FR", "NA"]
  types_simplified_name = types_name[:-1]
  
  num_type_classification = len(types_simplified_name)
  rxn_classification_num = [0]*(num_type_classification+1) #total number of classified cases for each type
  rxn_classification_num_PR = np.zeros((16, (num_type_classification+1))) #16 rows for 4prds*4rcts

  df_classification = pd.DataFrame(columns = COLUMN_NAME_df_classification)
  df_mol_stat = pd.DataFrame(columns = COLUMN_NAME_df_mol_stat)
  df_mol_stat_PR = {} #set of dataframes to save mol stat for each PR
  df_table_PR = pd.DataFrame(columns = ["R = 0", "R = 1", "R = 2", "R > 2"], \
                             index = ["P = 0", "P = 1", "P = 2", "P > 2"])
  df_table_PR_per_model = pd.DataFrame(columns = ["R = 0", "R = 1", "R = 2", "R > 2"], \
                             index = ["P = 0", "P = 1", "P = 2", "P > 2"])    

  for i in range(16):
    df_mol_stat_PR[i] = pd.DataFrame(columns = COLUMN_NAME_df_mol_stat)
  

  # for different input dataset
  iterator = simple_sbml.modelIterator(initial=initial_model_indx, final=final_model_indx,
  data_dir = data_dir, zip_filename = zip_filename)

  biomodel_non_count = 0
  for idx, item in enumerate(iterator):
    if item is None:
      file_num = initial_model_indx + idx
      print("File %d has an error." % (file_num))
    else:
      name = item.filename

      # Create an SBML model. We'll use the model
      # data/
      # try:
      #   path = os.path.join(cn.PROJECT_DIR, "data")
      # except:
      #   print("error")
      # path = os.path.join(path, name)
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
          #print("kinetics:", kinetics)

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

          #print("species_in_kinetic_law:", species_in_kinetic_law)
          #print("parameters_in_kinetic_law:", parameters_in_kinetic_law)

          #only for MM, MMcat and FR
          if len(reactant_list) != 0:
            ids_list += reactant_list # some rcts/prds also needs symbols definition
          if len(product_list) != 0:
            ids_list += product_list
          ids_list = list(dict.fromkeys(ids_list))

          # print("reactant_list:", reactant_list)
          #print("ids_list:", ids_list)

          #Define the keyword arguments
          kwargs = {"kinetics": kinetics, "kinetics_sim": kinetics_sim, \
            "reactant_list": reactant_list, "product_list": product_list, \
            "species_in_kinetic_law": species_in_kinetic_law, "parameters_in_kinetic_law": parameters_in_kinetic_law, \
            "ids_list": ids_list}

          classification_cp = [#needs to be in order
            reaction.kinetic_law.isZerothOrder(**kwargs),
            # reaction.kinetic_law.isPowerTerms(**kwargs),
            reaction.kinetic_law.isUNDR(**kwargs),
            reaction.kinetic_law.isUNMO(**kwargs),
            reaction.kinetic_law.isBIDR(**kwargs),
            reaction.kinetic_law.isBIMO(**kwargs),
            reaction.kinetic_law.isMM(**kwargs),
            reaction.kinetic_law.isMMcat(**kwargs),
            reaction.kinetic_law.isHill(**kwargs),
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
          
          if len(df_classification) == 0:
              df_classification = pd.DataFrame(classification_row_dct)
          else:
              df_classification = pd.concat([df_classification,\
                  pd.DataFrame(classification_row_dct)], ignore_index=True)


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

        if len(df_mol_stat) == 0:
            df_mol_stat = pd.DataFrame(mol_stat_row_dct)
        else:
            df_mol_stat = pd.concat([df_mol_stat,\
                pd.DataFrame(mol_stat_row_dct)], ignore_index=True)
      
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
            
            if len(df_mol_stat_PR[xy]) == 0:
                df_mol_stat_PR[xy] = pd.DataFrame(mol_stat_PR_row_dct)
            else:
                df_mol_stat_PR[xy] = pd.concat([df_mol_stat_PR[xy],\
                    pd.DataFrame(mol_stat_PR_row_dct)], ignore_index=True) 
  #for all the biomodels below here
  # This part is the same as the printed part in main section
  df_gen_stat = pd.DataFrame(columns = COLUMN_NAME_df_gen_stat)
  if(rxn_num != 0):
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

      if len(df_gen_stat) == 0:
          df_gen_stat = pd.DataFrame(gen_stat_row_dct)
      else:
          df_gen_stat = pd.concat([df_gen_stat,\
              pd.DataFrame(gen_stat_row_dct)], ignore_index=True) 
    
    df_gen_stat.at[0, RXN_NUM] = rxn_num
    df_gen_stat.at[0, BIOMOL_NUM] = len(df_mol_stat.index)


  #PR
  df_gen_stat_PR = pd.DataFrame(columns = COLUMN_NAME_df_gen_stat[0:-2])
  for xy in range(16):
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
        if len(df_gen_stat_PR) == 0:
            df_gen_stat_PR = pd.DataFrame(gen_stat_PR_row_dct)
        else:
            df_gen_stat_PR = pd.concat([df_gen_stat_PR,\
                pd.DataFrame(gen_stat_PR_row_dct)], ignore_index=True) 
    else: 
      gen_stat_PR_mul_row_dct = {k:[] for k in COLUMN_NAME_df_gen_stat[0:-2]}
      gen_stat_PR_mul_row_dct[CLASSIFICATIONS].extend(types_name)
      gen_stat_PR_mul_row_dct[PERCENTAGE].extend([0.]*(num_type_classification+1))
      gen_stat_PR_mul_row_dct[PERCENTAGE_PER_MODEL].extend([0.]*(num_type_classification+1))
      gen_stat_PR_mul_row_dct[PERCENTAGE_PER_MODEL_SDER].extend([0.]*(num_type_classification+1))

      if len(df_gen_stat_PR) == 0:
          df_gen_stat_PR = pd.DataFrame(gen_stat_PR_mul_row_dct)
      else:
          df_gen_stat_PR = pd.concat([df_gen_stat_PR,\
              pd.DataFrame(gen_stat_PR_mul_row_dct)], ignore_index=True) 

    df_table_PR.iloc[xy//4,xy%4] = rxn_num_PR[xy]
    if len(df_mol_stat_PR[xy]) != 0:
      df_table_PR_per_model.iloc[xy//4,xy%4] = rxn_num_PR[xy]/len(df_mol_stat_PR[xy].index)
    else:
      df_table_PR_per_model.iloc[xy//4,xy%4] = 0.

  return (df_classification, df_gen_stat, df_mol_stat, df_gen_stat_PR, biomodel_non_count, \
    df_table_PR, df_table_PR_per_model)


if __name__ == '__main__':
  start_time = time.time()

  initial_model_indx = 5
  final_model_indx = 6
  (df_classification, df_gen_stat, df_mol_stat, df_gen_stat_PR, biomodel_non_count, \
   df_table_PR, df_table_PR_per_model) = _dataSetStatistics(initial_model_indx = initial_model_indx, 
   final_model_indx = final_model_indx)
  rxn_num = len(df_classification)
  print(rxn_num)

  print("--- %s seconds ---" % (time.time() - start_time))
