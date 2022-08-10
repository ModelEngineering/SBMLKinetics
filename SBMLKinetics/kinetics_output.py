
# This script was written by Jin Xu and available on Github
# https://github.com/ModelEngineering/kinetics_validator
# This file includes all the functions to do the kinetics analysis.


from zipfile import ZIP_FILECOUNT_LIMIT
from isort import file
from SBMLKinetics import kinetics_classification
import sys

import numpy as np
import os

from sympy import *
from libsbml import * # access functions in SBML
import time

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib

from SBMLKinetics import types

# Column names
SBMLID = "SBMLid"
CLASSIFICATIONS = 'Classifications'
REACTION = 'Reaction'
KINETICLAW = 'kinetic law'
PERCENTAGE = 'Percentage'

class KineticAnalyzer:
  """
  Load Dataset of SBML files.

  Kinetics type (K type) including ten types "ZERO" (Zeroth order), "UNDR" 
  (Uni-directional mass action), "UNMO" (Uni-term with moderator), "BIDR" 
  (Bi-directional mass action), "BIMO" (Bi-terms with moderator), "MM" 
  (Michaelis-Menten kinetics without explicit enzyme), "MMCAT" 
  (Michaelis-Menten kinetics with explicit enzyme), "HILL" (Hill equations), 
  "FR" (Kinetic law in the format of fraction other than MM, MMCAT or HILL) and "NA" 
  (not classified kinetics). 

  Reaction type (R type) is quantitatively represented by the number of reactants 
  (r = 0, 1, 2, 3 (representing>2)) and products (p= 0, 1, 2, 3 (representing>2)).

  Args: 
      dataSet: str-"biomodels", "curated", "metabolic", "signalling", "homo_sapiens", "non_homo", 
      "cellular_organisms", "Mus_musculus", "Mammalia" or "Saccharomyces_cerevisiae";
      
      path: str-path to the file, with a format of ``D:\\path\\to``;
      
      model_indices: range-(initial_model_indx, final_model_indx)

  """

  def __init__(self, path = os.path.dirname(os.path.abspath(__file__)), 
    dataSet = "biomodels", model_indices = range(0,1000)):

    #In addition to dataSetName, allow users to inmport a zip of sbml files from a path 
    initial_model_indx = min(model_indices)
    final_model_indx = max(model_indices)

    if type(dataSet) == str and dataSet in ["biomodels", "curated", 
    "metabolic", "signalling", "homo_sapiens", "non_homo", 
    "cellular_organisms", "Mus_musculus", "Mammalia", "Saccharomyces_cerevisiae"]:
      zip_filename = dataSet + '.zip'
      try:
        self.tuple = kinetics_classification._dataSetStatistics(zip_filename = zip_filename, 
        initial_model_indx = initial_model_indx, final_model_indx = final_model_indx)
      except Exception as err:
          raise Exception (err)

    elif '.zip' in dataSet:
      try:
        self.tuple = kinetics_classification._dataSetStatistics(data_dir = path, zip_filename = dataSet, 
        initial_model_indx = initial_model_indx, final_model_indx = final_model_indx)
      except Exception as err:
          raise Exception (err)

    else:
      raise Exception("Not a valid dataset input.")

  ##Query Distributions    

  def getKTypeDistribution(self):
    """
    Get the kinetics type distribution.

    Returns:
        df_gen_stat_final: dataFrame-kinetics type distribution. 
        The column names are: "Classifications", "Percentage", "Percentage standard error", 
        "Percentage per model", "Percentage per model standard error".
        
        The column of "Classifications" covers the ten kinetic law types.
        
    """  
 
    (_, df_gen_stat, _, _, _, _, _) = self.tuple

    df_gen_stat_final = df_gen_stat[["Classifications", "Percentage", "Percentage per model", \
      "Percentage per model standard error"]]
    try:
      df_gen_stat_final.insert(2, "Percentage standard error", 0)
    except:
      pass

    return df_gen_stat_final


  def getKTypeDistributionPerRType(self, R_type):
    """
    Get the kinetics type distribution for the certein reaction type.

    Args: 
        R_type: an object with attributes rct_num and prd_num representing 
        the number of reactants and products.

        rct_num: int-0, 1, 2, 3 (representing > 2).
        
        prd_num: int-0, 1, 2, 3 (representing > 2).
        
    Returns:
        df_gen_stat_PR_final: dataFrame-the kinetics distribution for a reaction type.
        
        The column names are: "Classifications", "Percentage", "Percentage standard error", 
        "Percentage per model", "Percentage per model standard error".
        
        The column of "Classifications" covers the ten kinetic law types.
        
    """  
    rct_num = R_type.rct_num
    prd_num = R_type.prd_num

    (_, df_gen_stat, _, df_gen_stat_PR, _, _, _) = self.tuple

    df_gen_stat_plot = df_gen_stat[["Classifications", "Percentage", "Percentage per model", \
      "Percentage per model standard error"]]
    try:
      df_gen_stat_plot.insert(2, "Percentage standard error", 0)
    except:
      pass
    try:
      df_gen_stat_PR.insert(2, "Percentage standard error", 0)
    except:
      pass
    df_gen_stat_PR_plot = {}
    types = len(df_gen_stat_plot)
    if prd_num in [0,1,2,3] and rct_num in [0,1,2,3]:
      i = prd_num*4 + rct_num
      df_gen_stat_PR_plot[i] = pd.DataFrame(columns = df_gen_stat_PR.columns.tolist())
      df_temp = df_gen_stat_PR[types*i:types*(i+1)]   
      df_gen_stat_PR_plot[i] = pd.concat([df_gen_stat_PR_plot[i],df_temp], ignore_index=True)

      df_gen_stat_PR_final = df_gen_stat_PR_plot[i]
      return df_gen_stat_PR_final
    else:
      raise Exception("Not a valid reactant or product number.")

  def getRTypeDistribution(self):
    """
    Get the distribution of reaction involved for each reaction type.

    Returns:
        df_table_PR_final: dataFrame-Reaction type distribution. 
        
        The column names represent number of reactants: "R = 0", "R = 1", "R = 2", "R > 2".

        The row names represent number of products: "P = 0", "P = 1", "P = 2", "P > 2".
        
    """ 
    (_, _, _, _, _, df_table_PR, _) = self.tuple

    #generate the PR tables
    try:
      df_table_PR_final = df_table_PR.div(df_table_PR.sum().sum())

      return df_table_PR_final
    except Exception as e:
      raise Exception(e)

  def getRTypeDistributionPerModel(self):
    """
    Get the distribution of reaction involved for each reaction type per model.

    Returns:
        df_table_PR_per_model final: dataFrame-Reaction type distribution. 
        
        The column names represent number of reactants: "R = 0", "R = 1", "R = 2", "R > 2".

        The row names represent number of products: "P = 0", "P = 1", "P = 2", "P > 2".
        
    """ 
    (_, _, _, _, _, _, df_table_PR_per_model) = self.tuple

    #generate the PR tables
    try:
      df_table_PR_per_model_final = df_table_PR_per_model.div(df_table_PR_per_model.sum().sum())
      return df_table_PR_per_model_final
    except Exception as e:
      raise Exception(e)


  ##Query Elements
  def getTopKType(self):
    """
    Get the most frequent kinetics type. 

    Returns:
        kinetic_type_list: list of K_type. Sometimes there could be more than one
        top kinetics type to make the length of kinetic_type_list larger than one. 
        
        K_type: an object with an attribute K_type_str representing 
        the type of kinetics.

        K_type_str: str-"ZERO" (Zeroth order), "UNDR" 
        (Uni-directional mass action), "UNMO" (Uni-term with moderator), "BIDR" 
        (Bi-directional mass action), "BIMO" (Bi-terms with moderator), "MM" 
        (Michaelis-Menten kinetics without explicit enzyme), "MMCAT" 
        (Michaelis-Menten kinetics with explicit enzyme), "HILL" (Hill equations), 
        "FR" (Kinetic law in the format of fraction other than MM, MMCAT or HILL) and "NA" 
        (not classified kinetics). 
 
    """  

    df_temp = self.getKTypeDistribution()
    max_value = df_temp['Percentage'].max()
    idx_list = df_temp.index[df_temp['Percentage'] == max_value].tolist()
    kinetics_type_list =[] 
    for i in range(len(idx_list)):
        kinetics_type_list.append(types.K_type(df_temp.iloc[idx_list[i]]["Classifications"]))
      
    return kinetics_type_list
  
  def getKTypeProb(self, K_type):
    """
    Get the probability value of the certain kinetics type.

    Args:
        K_type: an object with an attribute K_type_str representing 
        the type of kinetic law.

        K_type_str: str-"ZERO" (Zeroth order), "UNDR" 
        (Uni-directional mass action), "UNMO" (Uni-term with moderator), "BIDR" 
        (Bi-directional mass action), "BIMO" (Bi-terms with moderator), "MM" 
        (Michaelis-Menten kinetics without explicit enzyme), "MMCAT" 
        (Michaelis-Menten kinetics with explicit enzyme), "HILL" (Hill equations), 
        "FR" (Kinetic law in the format of fraction other than MM, MMCAT or HILL) and "NA" 
        (not classified kinetics). 

    Returns:
        kinetics_value: float-the probability of the certain kinetics type.
    """  
    K_type_str = K_type.K_type_str
    df_temp = self.getKTypeDistribution()
    idx_list = df_temp.index[df_temp['Classifications'] == K_type_str].tolist()
    if len(idx_list) == 0:
      raise Exception("Please enter a valid kinetic type.")
    else:
      kinetics_value_list =[] 
      for i in range(len(idx_list)):
          kinetics_value_list.append(df_temp.iloc[idx_list[i]]["Percentage"])
      kinetics_value = kinetics_value_list[0] 
      
      return kinetics_value

  def getTopKTypePerRType(self, R_type):

    """
    Get the most frequent kinetics type from a certain reaction type. 

    Args: 
        R_type: an object with attributes rct_num and prd_num representing 
        the number of reactants and products.

        rct_num: int-0, 1, 2, 3 (representing > 2).
        
        prd_num: int-0, 1, 2, 3 (representing > 2).

    Returns:
        kinetic_type_list: list of k_type. Sometimes there could be more than one
        top kinetics type to make the length of kinetic_type_list larger than one. 
         
        K_type: an object with an attribute K_type_str representing 
        the type of kinetic law.

        K_type_str: str-"ZERO" (Zeroth order), "UNDR" 
        (Uni-directional mass action), "UNMO" (Uni-term with moderator), "BIDR" 
        (Bi-directional mass action), "BIMO" (Bi-terms with moderator), "MM" 
        (Michaelis-Menten kinetics without explicit enzyme), "MMCAT" 
        (Michaelis-Menten kinetics with explicit enzyme), "HILL" (Hill equations), 
        "FR" (Kinetic law in the format of fraction other than MM, MMCAT or HILL) and "NA" 
        (not classified kinetics). 
    """  
    rct_num = R_type.rct_num
    prd_num = R_type.prd_num

    if prd_num in [0,1,2,3] and rct_num in [0,1,2,3]:
      df_temp = self.getKTypeDistributionPerRType(R_type = R_type)

      max_value = df_temp['Percentage'].max()
      idx_list = df_temp.index[df_temp['Percentage'] == max_value].tolist()
      kinetics_type_list =[] 
      for i in range(len(idx_list)):
          kinetics_type_list.append(types.K_type(df_temp.iloc[idx_list[i]]["Classifications"]))
      
      return kinetics_type_list

    else:
      raise Exception("Not a valid reactant or product number.")

  def getKTypeProbPerRType(self, R_type, K_type):
    """
    Get the probability value of the certain kinetics type from a certain reaction type.

    Args:
        R_type: an object with attributes rct_num and prd_num representing 
        the number of reactants and products.

        rct_num: int-0, 1, 2, 3 (representing > 2).
        
        prd_num: int-0, 1, 2, 3 (representing > 2).

        K_type: an object with an attribute K_type_str representing the type of kinetics.

        K_type_str: str-"ZERO" (Zeroth order), "UNDR" 
        (Uni-directional mass action), "UNMO" (Uni-term with moderator), "BIDR" 
        (Bi-directional mass action), "BIMO" (Bi-terms with moderator), "MM" 
        (Michaelis-Menten kinetics without explicit enzyme), "MMCAT" 
        (Michaelis-Menten kinetics with explicit enzyme), "HILL" (Hill equations), 
        "FR" (Kinetic law in the format of fraction other than MM, MMCAT or HILL) and "NA" 
        (not classified kinetics). 

    Returns:
        kinetics_value: float-the probability of the certain kinetics type.
    """  
    K_type_str = K_type.K_type_str

    df_temp = self.getKTypeDistributionPerRType(R_type = R_type)

    idx_list = df_temp.index[df_temp['Classifications'] == K_type_str].tolist()
    if len(idx_list) == 0:
      raise Exception("Please enter a valid kinetic type.")
    else:
      kinetics_value_list =[] 
      for i in range(len(idx_list)):
          kinetics_value_list.append(df_temp.iloc[idx_list[i]]["Percentage"])
      kinetics_value = kinetics_value_list[0] 
      
      return kinetics_value

  def getTopRType(self):

    """
    Get the most frequent reaction type (with the most number of reactions involved in
    the certain reaction type). 

    returns: 
        rct_prd_num_list: list of R_type info with the most frequent reaction 
        type. Sometimes there could be more than one top kinetics type to make the length 
        of rct_prd_num_list larger than one. 
        
        R_type: an object with attributes rct_num and prd_num representing 
        the number of reactants and products.

        rct_num: int-0, 1, 2, 3 (representing > 2).
        
        prd_num: int-0, 1, 2, 3 (representing > 2).

    """  

    df_temp = self.getRTypeDistributionPerModel()

    max_value = 0
    rct_prd_num_list = []

    for col in df_temp.columns:
      max = df_temp[col].max()
      idx_list = df_temp.index[df_temp[col] == max].tolist()
      if max > max_value:
        max_value = max

    for col in df_temp.columns:
      idx_list = df_temp.index[df_temp[col] == max_value].tolist()
      for i in range(len(idx_list)):
        if '=' in col:
          rct_num = int(col[-1])
        else:
          rct_num = 3
        if '=' in idx_list[i]:
          prd_num = int(idx_list[i][-1])
        else:
          prd_num = 3 
        
        rct_prd_num_list.append(types.R_type(rct_num, prd_num))

    return rct_prd_num_list

  def getRTypeProb(self, R_type):

    """
    Args:
        R_type: an object with attributes rct_num and prd_num representing 
        the number of reactants and products.

        rct_num: int-0, 1, 2, 3 (representing > 2).
        
        prd_num: int-0, 1, 2, 3 (representing > 2).

    Returns:
        R_value: float-the probability of the certain reaction type.
    """  
    rct_num = R_type.rct_num
    prd_num = R_type.prd_num

    if prd_num in [0,1,2,3] and rct_num in [0,1,2,3]:
      df_temp = self.getRTypeDistribution()
      value = df_temp.iat[prd_num, rct_num]
    else:
        raise Exception("Not a valid reactant or product number.")

    return value

  def getNumBiomodelsAnalyzed(self):
    """
    Get the number of biomodels analyzed.
    
    Returns:
        biomodels_num: int-number of biomodels.
    """  
    (_, _, df_mol_stat, _, _, _, _) = self.tuple
    biomodels_num = len(df_mol_stat)
    
    return biomodels_num
  

  def getNumRxnsAnalyzed(self):
    """
    Get the number of reactions analyzed.
    
    Returns:
        rxn_num: int-number of reactions.
    """  
    (df_classification, _, _, _, _, _, _) = self.tuple
    rxn_num = len(df_classification)
    
    return rxn_num


  ##Plots 
  def plotKTypeDistribution(self, path = "", fileName = 'KTypeDistribution.pdf'):
    """
    Plot the kinetics type distribution and save it as a pdf file.

    Args: 
        path: str-path to the file, with a format like ``D:/path/to/`` (or ``D:\\\path\\\ to\\\``)

        fileName: str-file name with which the pdf file save to, ending with ".pdf".
        
    """  
    if str(fileName).lower()[-4:] == ".pdf" and fileName[:-4] != "":
      df_gen_stat_plot = self.getKTypeDistribution()
      path_fileName = path + fileName[:-4] + ".pdf"

      yerr = df_gen_stat_plot[["Percentage standard error", \
        "Percentage per model standard error"]].to_numpy().T
      ax = df_gen_stat_plot.plot(kind="bar",x="Classifications", y=["Percentage","Percentage per model"],\
        yerr=yerr, fontsize = 8)
      ax.set_ylim(0.,1.)
      ax.get_yaxis().set_major_formatter(
        matplotlib.ticker.FuncFormatter(lambda y, p: str("{:.2%}".format(y))))
      for p in ax.patches:
        ax.annotate(str("{:.2%}".format(p.get_height())), (p.get_x() * 1.005, p.get_height() * 1.005), fontsize = 4)
      #plt.show()
      fig = ax.get_figure()
      fig.savefig(path_fileName)
    else:
      raise Exception("Please enter a valid pdf file name.")


  def plotKTypeDistributionPerRType(self, R_type, path = "", fileName = "KTypeDistributionPerRType.pdf"):
    """
    Plot the kinetics type distribution for a certain reaction type and save it as
    a pdf file.

    Args:  
        R_type: an object with attributes rct_num and prd_num representing 
        the number of reactants and products.

        rct_num: int-0, 1, 2, 3 (representing > 2).
        
        prd_num: int-0, 1, 2, 3 (representing > 2).

        path: str-path to the file, with a format like ``D:/path/to/`` (or ``D:\\\path\\\ to\\\``)

        fileName: str-file name with which the pdf file save to, ending with ".pdf".
    """  

    rct_num = R_type.rct_num
    prd_num = R_type.prd_num

    if str(fileName).lower()[-4:] == ".pdf" and fileName[:-4] != "":
      path_fileName = path + fileName[:-4] + ".pdf"
      (_, _, _, _, _, \
      df_table_PR, df_table_PR_per_model) = self.tuple

      #generate the PR two tables
      try:
        df_table_PR_plot = df_table_PR.div(df_table_PR.sum().sum())
        df_table_PR_per_model_plot = df_table_PR_per_model.div(df_table_PR_per_model.sum().sum())
      except Exception as e:
        raise Exception(e)

      df_gen_stat_PR_plot = {}

      if prd_num in [0,1,2,3] and rct_num in [0,1,2,3]:
        i = prd_num*4 + rct_num
        df_gen_stat_PR_plot[i] = self.getKTypeDistributionPerR(R_type = R_type)
        yerr = df_gen_stat_PR_plot[i][["Percentage standard error", \
          "Percentage per model standard error"]].to_numpy().T
        ax = df_gen_stat_PR_plot[i].plot(kind="bar",x="Classifications", y=["Percentage","Percentage per model"],\
          yerr=yerr, legend = None, fontsize = 8)
        ax.set_ylim(0.,1.)
        ax.get_yaxis().set_major_formatter(
          matplotlib.ticker.FuncFormatter(lambda y, p: str("{:.2%}".format(y))))
        for p in ax.patches:
          ax.annotate(str("{:.2%}".format(p.get_height())), (p.get_x() * 1.005, p.get_height() * 1.005), fontsize = 4)
        #plt.show()
        ax.get_yaxis().set_major_formatter(
        matplotlib.ticker.FuncFormatter(lambda y, p: str("{:.2%}".format(y))))
        ax.annotate('%s'%"{:.2%}".format(df_table_PR_plot.iat[i//4, i%4]), xy=(0, .9), color = 'dodgerblue')
        if i//4 == 3 and i % 4 != 3:
          ax.annotate('P > %d, R = %d'%(2, i%4), xy=(3, .9))
        elif i//4 != 3 and i % 4 == 3:
          ax.annotate('P = %d, R > %d'%(i//4, 2), xy=(3, .9))
        elif i//4 == 3 and i % 4 == 3:
          ax.annotate('P > %d, R > %d'%(2, 2), xy=(3, .9))
        else:
          ax.annotate('P = %d, R = %d'%(i//4, i%4), xy=(3, .9))
        ax.annotate('%s'%"{:.2%}".format(df_table_PR_per_model_plot.iat[i//4, i%4]), xy=(7., .9), color = 'darkorange')
        fig = ax.get_figure()
        handles, labels = fig.axes[-1].get_legend_handles_labels()
        fig.legend(handles, labels, loc='upper center')
        fig.savefig(path_fileName)
      else:
        raise Exception("Not a valid reactant or product number.")
    else:
      raise Exception("Please enter a valid pdf file name.")


  def plotKTypeDistributionVsRType(self, path = "", fileName = 'KTypeDistributionVsRType.pdf'):
    """
    Plot the kinetics type distribution vs each reaction type and save it as a pdf 
    file.
  
    Args: 
        path: str-path to the file, with a format like ``D:/path/to/`` (or ``D:\\\path\\\ to\\\``)

        fileName: str-file name with which the pdf file save to, ending with '.pdf'.

    """  
    
    if str(fileName).lower()[-4:] == ".pdf" and fileName[:-4] != "":
      path_fileName = path + fileName[:-4] + ".pdf"
      (_, _, _, _, _, \
      df_table_PR, df_table_PR_per_model) = self.tuple

      #generate the PR two tables
      try:
        df_table_PR_plot = df_table_PR.div(df_table_PR.sum().sum())
        df_table_PR_per_model_plot = df_table_PR_per_model.div(df_table_PR_per_model.sum().sum())
      except Exception as e:
        raise Exception(e)

      df_gen_stat_PR_plot = {}

      fig = plt.figure(figsize = (16,16))
      axes = fig.subplots(nrows=4, ncols=4)
      for i in range(16):
        df_gen_stat_PR_plot[i] = \
          self.getKTypeDistributionPerRType(types.R_type(rct_num=i//4, prd_num=i%4))
        yerr = df_gen_stat_PR_plot[i][["Percentage standard error", \
          "Percentage per model standard error"]].to_numpy().T
        df_gen_stat_PR_plot[i].plot(ax = axes[i//4,i%4] , kind="bar", 
            x="Classifications", y=["Percentage","Percentage per model"],\
            yerr=yerr, legend = None, fontsize = 6)
        axes[i//4, i%4].get_yaxis().set_major_formatter(
        matplotlib.ticker.FuncFormatter(lambda y, p: str("{:.2%}".format(y))))
        axes[i//4, i%4].annotate('%s'%"{:.2%}".format(df_table_PR_plot.iat[i//4, i%4]), xy=(0, .9), color = 'dodgerblue')
        if i//4 == 3 and i % 4 != 3:
          axes[i//4, i%4].annotate('P > %d, R = %d'%(2, i%4), xy=(3, .9))
        elif i//4 != 3 and i % 4 == 3:
          axes[i//4, i%4].annotate('P = %d, R > %d'%(i//4, 2), xy=(3, .9))
        elif i//4 == 3 and i % 4 == 3:
          axes[i//4, i%4].annotate('P > %d, R > %d'%(2, 2), xy=(3, .9))
        else:
          axes[i//4, i%4].annotate('P = %d, R = %d'%(i//4, i%4), xy=(3, .9))
        axes[i//4, i%4].annotate('%s'%"{:.2%}".format(df_table_PR_per_model_plot.iat[i//4, i%4]), xy=(7., .9), color = 'darkorange')
        #if i//4 != 3:
        if i != 12:
          axes[i//4, i%4].get_xaxis().set_visible(False)
        axes[i//4,i%4].set_ylim([0, 1])  
      handles, labels = fig.axes[-1].get_legend_handles_labels()
      fig.legend(handles, labels, loc='upper center')
      #plt.show()
      fig.savefig(path_fileName)
    else:
      raise Exception("Please enter a valid pdf file name.")


  def plotRTypeDistribution(self, path = "", fileName = 'RTypeDistribution.pdf'):
    """
    Plot the distribution of reaction involved for reaction type and save it 
    as a pdf file.
  
    Args: 
        path: str-path to the file, with a format like ``D:/path/to/`` (or ``D:\\\path\\\ to\\\``)

        fileName: str-file name with which the pdf file save to, ending with '.pdf'.
    """  
    if str(fileName).lower()[-4:] == ".pdf" and fileName[:-4] != "":
      path_fileName = path + fileName[:-4] + ".pdf"

      df_table_PR_plot = self.getRTypeDistribution()

      # Defining figure size for the output plot 
      fig, bx = plt.subplots(figsize = (12, 7))
      idx = df_table_PR_plot.index.tolist()
      cols = df_table_PR_plot.columns.tolist()
      df_table_PR_plot = pd.DataFrame(df_table_PR_plot.values.tolist(),
                      columns = cols, index = idx)   
      # Displaying dataframe as an heatmap 
      # with diverging colourmap as RdYlGn
      fmt = lambda x, pos: '{:.2%}'.format(x)
      sns.heatmap(df_table_PR_plot, cmap ='RdYlGn', linewidths = 0.30, annot = True, fmt='.2%')
      cbar = bx.collections[0].colorbar
      cbar.set_ticklabels(['0', '10%', '20%', '30%', '40%', '50%', '60%', '70%', '80%', '90%'])
      # Displaying the figure
      #plt.show()  
      plt.savefig(path_fileName, dpi=350)
    else:
      raise Exception("Please enter a valid pdf file name.")


  def plotRTypeDistributionPerModel(self, path = "", fileName = 'RTypeDistributionPerModel.pdf'):
    """
    Plot the distribution of reactions involved for each reaction type per model and 
    save it as a pdf file.
      
    Args: 
        path: str-path to the file, with a format like ``D:/path/to/`` (or ``D:\\\path\\\ to\\\``)

        fileName: str-file name with which the pdf file save to, ending with '.pdf'.
    """  
    if str(fileName).lower()[-4:] == ".pdf" and fileName[:-4] != "":
      path_fileName = path + fileName[:-4] + ".pdf"

      df_table_PR_per_model_plot = self.getRTypeDistributionPerModel()

      fig, cx = plt.subplots(figsize = (12, 7))
      idx = df_table_PR_per_model_plot.index.tolist()
      cols = df_table_PR_per_model_plot.columns.tolist()
      df_table_PR_per_model_plot = pd.DataFrame(df_table_PR_per_model_plot.values.tolist(),
                      columns = cols, index = idx)   
      # Displaying dataframe as an heatmap 
      # with diverging colourmap as RdYlGn
      fmt = lambda x, pos: '{:.2%}'.format(x)
      sns.heatmap(df_table_PR_per_model_plot, cmap ='RdYlGn', linewidths = 0.30, annot = True, fmt='.2%')
      cbar = cx.collections[0].colorbar
      cbar.set_ticklabels(['0', '5%', '10%', '15%', '20%', '25%', '30%', '35%', \
        '40%', '45%', '50%', '55%', '60%', '65%', '70%', '75%', '80%', '85%', '90%', '95%'])
      # Displaying the figure
      #plt.show()  
      plt.savefig(path_fileName, dpi=350)
    else:
      raise Exception("Please enter a valid pdf file name.")
  
  def _tableKTypeDistribution(self, path = "", fileName = "KTypeDistribution.xlsx"):
    """
    Save the kinetics type distribution to an excel file.

    Args: 
        path: str-path to the file, with a format like ``D:/path/to/`` (or ``D:\\\path\\\ to\\\``)

        fileName: str-file name with which the excel file save to, ending with ".xlsx".
        
    """  

    if str(fileName).lower()[-5:] == ".xlsx" and fileName[:-5] != "":
      df_gen_stat_final = self.getKTypeDistribution()
      # Create a Pandas Excel writer using XlsxWriter as the engine.
      path_fileName = path + fileName[:-5] + ".xlsx"
      writer = pd.ExcelWriter(path_fileName, engine='xlsxwriter')
      df_gen_stat_final.to_excel(writer, sheet_name='general_statistics')
      # Close the Pandas Excel writer and output the Excel file.
      writer.save()
    else:
      raise Exception("Please enter a valid excel file name.")

  def _tableKTypeDistributionPerRType(self, R_type, path = "", fileName = "KTypeDistributionPerRType.xlsx"):
    """
    Save the kinetics type distribution for a certain reaction type to an excel file.

    Args: 
        R_Type: an object with attributes rct_num and prd_num representing 
        the number of reactants and products.

        rct_num: int-0, 1, 2, 3 (representing > 2).
        
        prd_num: int-0, 1, 2, 3 (representing > 2).

        path: str-path to the file, with a format like ``D:/path/to/`` (or ``D:\\\path\\\ to\\\``)
        
        fileName: str-file name with which the excel file save to, ending with ".xlsx".
        
    """

    if str(fileName).lower()[-5:] == ".xlsx" and fileName[:-5] != "":
      df_gen_stat_final = self.getKTypeDistributionPerRType(R_type = R_type)
      # Create a Pandas Excel writer using XlsxWriter as the engine.
      path_fileName = path + fileName[:-5] + ".xlsx"
      writer = pd.ExcelWriter(path_fileName, engine='xlsxwriter')
      df_gen_stat_final.to_excel(writer, sheet_name='general_statistics')
      # Close the Pandas Excel writer and output the Excel file.
      writer.save()
    else:
      raise Exception("Please enter a valid excel file name.")

  def _tableRTypeDistribution(self, path = "", fileName = "RTypeDistribution.xlsx"):
    """
    Save the distribution of reaction involved for each reaction type to an excel file.

    Args: 
        path: str-path to the file, with a format like ``D:/path/to/`` (or ``D:\\\path\\\ to\\\``)

        fileName: str-file name with which the excel file save to, ending with ".xlsx".
        
    """  

    if str(fileName).lower()[-5:] == ".xlsx" and fileName[:-5] != "":
      df_gen_stat_final = self.getRTypeDistribution()
      # Create a Pandas Excel writer using XlsxWriter as the engine.
      path_fileName = path + fileName[:-5] + ".xlsx"
      writer = pd.ExcelWriter(path_fileName, engine='xlsxwriter')
      df_gen_stat_final.to_excel(writer, sheet_name='general_statistics')
      # Close the Pandas Excel writer and output the Excel file.
      writer.save()
    else:
      raise Exception("Please enter a valid excel file name.") 

  def _tableRTypeDistributionPerModel(self, path = "", fileName = "RTypeDistributionPerModel.xlsx"):
    """
    Save the distribution of reaction involved for each reaction type per model to an 
    excel file.

    Args: 
        path: str-path to the file, with a format like ``D:/path/to/`` (or ``D:\\\path\\\ to\\\``)

        fileName: str-file name with which the excel file save to, ending with ".xlsx".
        
    """  

    if str(fileName).lower()[-5:] == ".xlsx" and fileName[:-5] != "":
      df_gen_stat_final = self.getRTypeDistributionPerModel()
      # Create a Pandas Excel writer using XlsxWriter as the engine.
      path_fileName = path + fileName[:-5] + ".xlsx"
      writer = pd.ExcelWriter(path_fileName, engine='xlsxwriter')
      df_gen_stat_final.to_excel(writer, sheet_name='general_statistics')
      # Close the Pandas Excel writer and output the Excel file.
      writer.save()
    else:
      raise Exception("Please enter a valid excel file name.") 

  def _printBriefStatOfKTypeDistribution(self):
    """
    Print the brief statistics for the kinetics type distribution.
    """  
    (df_classification, df_gen_stat, _, _, biomodel_non_count, _, _) = self.tuple
    rxn_num = len(df_classification)
    if(rxn_num != 0):
      #print("\n\n")
      print("A brief statistics for the classification of reactions:")
      print("Reaction number:", rxn_num)
      for i in range(len(df_gen_stat)):
        print(df_gen_stat.iloc[i][CLASSIFICATIONS] + ":" + str(df_gen_stat.iloc[i][PERCENTAGE]))
    else:
      print("There are no reactions.")
    print("number of biomodels with some reactions not classified:", biomodel_non_count)

  def _printKTypePerRxn(self):
    """
    Print the kinetics type for each reaction.
    """  
    (df_classification, _, _, _, _, _, _) = self.tuple
    SBML_id_list = []
    for i in range(len(df_classification)):
      SBML_id = df_classification.iloc[i][SBMLID]
      if SBML_id not in SBML_id_list:
        print(SBML_id)
        SBML_id_list.append(SBML_id)
      print(df_classification.iloc[i][REACTION] + ";" + df_classification.iloc[i][KINETICLAW])
      print(df_classification.iloc[i][CLASSIFICATIONS])

  def _saveAllStatisticsInfoToExcel(self, fileName):
    """
    Save all the statistics information of kinetics classification into an excel file.
      
    Args: 
        fileName: str-file name with which the excel file save to.

    """  
    (df_classification, df_gen_stat, df_mol_stat, df_gen_stat_PR, _, \
    df_table_PR, df_table_PR_per_model) = self.tuple

    #generate the PR two tables
    try:
      df_table_PR_plot = df_table_PR.div(df_table_PR.sum().sum())
      df_table_PR_per_model_plot = df_table_PR_per_model.div(df_table_PR_per_model.sum().sum())
    except Exception as e:
      raise Exception(e)

    # Create a Pandas Excel writer using XlsxWriter as the engine.
    writer = pd.ExcelWriter(fileName, engine='xlsxwriter')
    # Write each dataframe to a different worksheet.
    df_classification.to_excel(writer, sheet_name='classification')
    df_gen_stat.to_excel(writer, sheet_name='general_statistics')
    df_mol_stat.to_excel(writer, sheet_name='statistics_per_model')
    df_gen_stat_PR.to_excel(writer, sheet_name='general_statistics_PR')
    df_table_PR_plot.to_excel(writer, sheet_name = 'table_PR')
    df_table_PR_per_model_plot.to_excel(writer, sheet_name = 'table_PR_per_model')
    # Close the Pandas Excel writer and output the Excel file.
    writer.save()

if __name__ == '__main__':
  start_time = time.time()

  initial_model_indx = 5
  final_model_indx = 6

  model_indices = range(initial_model_indx, final_model_indx+1)
  analyzer = KineticAnalyzer(path = 'D:/summer-2020/Jo/kinetics_validator/data',
  dataSet = "biomodels.zip", model_indices=model_indices) 

  # #Query Distributions
  # print(analyzer.getKTypeDistribution()) 
  # print(analyzer.getKTypeDistributionPerRType(R_type = types.R_type(1,1)))
  # print(analyzer.getRTypeDistribution())
  # print(analyzer.getRTypeDistributionPerModel())

  
  #Query Elements
  # print(analyzer.getTopKType()[0].K_type_str)
  # print(analyzer.getKTypeProb(K_type = types.K_type("NA")))
  # print(analyzer.getTopKTypePerRType(R_type = types.R_type(1,1))[0].K_type_str)
  # print(analyzer.getKTypeProbPerRType(R_type = types.R_type(1,1), K_type=types.K_type("NA")))
  # print(analyzer.getTopRType()[0].rct_num)
  # print(analyzer.getRTypeProb(R_type = types.R_type(1,1)))
  # print(analyzer.getNumBiomodelsAnalyzed())
  # print(analyzer.getNumRxnsAnalyzed())


  # #Plots #tests are not necessary
  # analyzer.plotKTypeDistribution(path = 'D:/summer-2020/Jo/')
  # analyzer.plotKTypeDistributionPerRType(R_type = types.R_type(1,1))
  # analyzer.plotKTypeDistributionVsRType()
  # analyzer.plotRTypeDistribution()
  # analyzer.plotRTypeDistributionPerModel()

  # #Tables
  # analyzer._tableKTypeDistribution(path = 'D:/summer-2020/Jo/', fileName = "Kinetics.xlsx")
  # analyzer._tableKTypeDistribution()
  # analyzer._tableKTypeDistributionPerRType(R_type = types.R_type(1,1))
  # analyzer._tableRTypeDistribution()
  # analyzer._tableRTypeDistributionPerModel()


  # analyzer._printBriefStatOfKTypeDistribution()
  # analyzer._printKTypePerRxn()
  # analyzer._saveAllStatisticsInfoToExcel(fileName='statistics_result.xlsx')
 

  print("--- %s seconds ---" % (time.time() - start_time))
