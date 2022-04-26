
#This script is to do kinetic classification.
#Make sure that you have setup your PYTHONPATH environment
#variable as described in the github repository.


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

# Column names
SBMLID = "SBMLid"
CLASSIFICATIONS = 'Classifications'
REACTION = 'Reaction'
KINETICLAW = 'kinetic law'
PERCENTAGE = 'Percentage'

class KineticAnalyzer:
  """
  Load Dataset of SBML files.

  Args: 
      dataSet: str-"biomodels", "curated", "metabolic", "signalling", "homo_sapiens", "non_homo", 
      "cellular_organisms", "Mus_musculus", "Mammalia", "Saccharomyces_cerevisiae";
      
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

  def getKineticLawDistribution(self, path = "", fileName = ""):
    """
    Get the kinetic law distribution (and save the dataframe into an excel file).

    Args: 
        path: str-path to the file, with a format like ``D:/path/to/`` (or ``D:\\\path\\\ to\\\``)

        fileName: str-file name with which the excel file save to, "" (do not save to excel file).

    Returns:
        df_gen_stat_final: dataFrame-kinetic law distribution. 
        The column names are: "Classifications", "Percentage", "Percentage standard error", 
        "Percentage per model", "Percentage per model standard error".
        
        In the column of "Classifications", there are "ZERO", "UNDR", "UNMO", "BIDR", "BIMO", 
        "MM", "MMCAT", "HILL", "FR" and "NA" in detail. 
        
        "ZERO" means "Zeroth order", "UNDR" means "Uni-directional mass action", "UNMO" means
        "Uni-term with moderator", "BIDR" means "Bi-directional mass action", "BIMO" means "Bi-
        terms with moderator", "MM" means "Michaelis-Menten kinetics", "MMCAT" means "Michaelis-
        Menten kinetics", "HILL" means "Hill equations", "FR" means kinetics in the format of 
        fraction other than MM, MMCAT and HILL, "NA" means not classified kinetics.

    """  
 
    (_, df_gen_stat, _, _, _, _, _) = self.tuple

    df_gen_stat_final = df_gen_stat[["Classifications", "Percentage", "Percentage per model", \
      "Percentage per model standard error"]]
    try:
      df_gen_stat_final.insert(2, "Percentage standard error", 0)
    except:
      pass
    if fileName != "":
      # Create a Pandas Excel writer using XlsxWriter as the engine.
      path_fileName = path + fileName
      writer = pd.ExcelWriter(path_fileName, engine='xlsxwriter')
      df_gen_stat_final.to_excel(writer, sheet_name='general_statistics')
      # Close the Pandas Excel writer and output the Excel file.
      writer.save()

    return df_gen_stat_final

  def TopFrequentKineticLawType(self):
    """
    Return the most frequent kinetic law type on average in the loaded SBML dataset . 

    Returns:
        kinetics_type_list: list pf kinetics_type.
        
        kinetics_type: str-kinetic law type.  
    """  
 
    (_, df_gen_stat, _, _, _, _, _) = self.tuple

    df_gen_stat_plot = df_gen_stat[["Classifications", "Percentage", "Percentage per model", \
      "Percentage per model standard error"]]
    df_temp = df_gen_stat_plot

    # try:
    #   kinetics_type_list = []
    #   max_idx = df_temp['Percentage'].idxmax()
    #   kinetics_type = df_temp['Classifications'][max_idx]
    #   kinetics_type_list.append(kinetics_type)
    # except:
    max_value = df_temp['Percentage'].max()
    idx_list = df_temp.index[df_temp['Percentage'] == max_value].tolist()
    kinetics_type_list =[] 
    for i in range(len(idx_list)):
        kinetics_type_list.append(df_temp.iloc[idx_list[i]]["Classifications"])
      
    return kinetics_type_list


  def plotKineticLawDistribution(self, fileName = 'KineticLawDistribution.pdf'):
    """
    Plot the kinetic law distribution as save it as a pdf file.

    Args: 
        fileName: str-file name with which the pdf file save to.

    """  
    (_, df_gen_stat, _, _, _, _, _) = self.tuple

    df_gen_stat_plot = df_gen_stat[["Classifications", "Percentage", "Percentage per model", \
      "Percentage per model standard error"]]
    try:
      df_gen_stat_plot.insert(2, "Percentage standard error", 0)
    except:
      pass
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
    fig.savefig(fileName)


  def getKineticLawDistributionPerMassTransfer(self, rct_num, prd_num, path = "", fileName = ""):
    """
    Get the kinetic law distribution for the certein mass transfer 
    (and save the dataframe into an excel file).

    Args: 
        rct_num: int-0, 1, 2, 3 (representing > 2).
        
        prd_num: int-0, 1, 2, 3 (representing > 2).
        
        path: str-path to the file, with a format like ``D:/path/to/`` (or ``D:\\\path\\\ to\\\``)
        
        fileName: str-file name with which the excel file save to, "" (do not save to excel file).

    Returns:
        df_gen_stat_PR_final: dataFrame-the kinetic law distribution for a certain mass trasfer.
        
        The column names are: "Classifications", "Percentage", "Percentage standard error", 
        "Percentage per model", "Percentage per model standard error".
        
        In the column of "Classifications", there are "ZERO", "UNDR", "UNMO", "BIDR", "BIMO", 
        "MM", "MMCAT", "HILL", "FR" and "NA" in detail. 
        
        "ZERO" means "Zeroth order", "UNDR" means "Uni-directional mass action", "UNMO" means
        "Uni-term with moderator", "BIDR" means "Bi-directional mass action", "BIMO" means "Bi-
        terms with moderator", "MM" means "Michaelis-Menten kinetics", "MMCAT" means "Michaelis-
        Menten kinetics", "HILL" means "Hill equations", "FR" means kinetics in the format of 
        fraction other than MM, MMCAT and HILL, "NA" means not classified kinetics.
    """  
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

      if fileName != "":
        # Create a Pandas Excel writer using XlsxWriter as the engine.
        path_fileName = path + fileName
        writer = pd.ExcelWriter(path_fileName, engine='xlsxwriter')
        df_gen_stat_PR_plot[i].to_excel(writer, sheet_name='general_statistics_PR')
        # Close the Pandas Excel writer and output the Excel file.
        writer.save()

      df_gen_stat_PR_final = df_gen_stat_PR_plot[i]
      return df_gen_stat_PR_final
    else:
      raise Exception("Not a valid reactant or product number.")

  def TopFrequentKineticLawTypePerMassTransfer(self, rct_num, prd_num):

    """
    Return the most frequent kinetic law type on average for a certain mass transfer 
    in the loaded SBML dataset . 

    Args: 
        rct_num: int-0, 1, 2, 3 (representing > 2).
        
        prd_num: int-0, 1, 2, 3 (representing > 2).

    Returns:
        kinetics_type_list: list pf kinetics_type.
        
        kinetics_type: str-kinetic law type.
    """  
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

      df_temp = df_gen_stat_PR_plot[i]
      # try:
      #   kinetics_type_list = []
      #   max_idx = df_temp['Percentage'].idxmax()
      #   kinetics_type = df_temp['Classifications'][max_idx]
      #   kinetics_type_list.append(kinetics_type)
      # except:
      max_value = df_temp['Percentage'].max()
      idx_list = df_temp.index[df_temp['Percentage'] == max_value].tolist()
      kinetics_type_list =[] 
      for i in range(len(idx_list)):
          kinetics_type_list.append(df_temp.iloc[idx_list[i]]["Classifications"])
      
      return kinetics_type_list

    else:
      raise Exception("Not a valid reactant or product number.")

  def plotKineticLawDistributionPerMassTransfer(self, rct_num, prd_num, 
  fileName = "KineticLawDistributionPerMassTransfer.pdf"):
    """
    Plot the kinetic law distribution for the certain mass transfer.

    Args: 
        rct_num: int - 0, 1, 2, 3 (representing > 2)
        
        prd_num: int - 0, 1, 2, 3 (representing > 2)
        
        fileName: str-file name with which the pdf file save to. 

    """  
    (_, df_gen_stat, _, df_gen_stat_PR, _, \
    df_table_PR, df_table_PR_per_model) = self.tuple

    #generate the PR two tables
    try:
      df_table_PR_plot = df_table_PR.div(df_table_PR.sum().sum())
      df_table_PR_per_model_plot = df_table_PR_per_model.div(df_table_PR_per_model.sum().sum())
    except Exception as e:
      raise Exception(e)

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
      fig.savefig(fileName)
    else:
      raise Exception("Not a valid reactant or product number.")


  def plotKineticLawDistributionVsMassTransfer(self, fileName = 'KineticLawDistributionVsMassTransfer.pdf'):
    """
    Plot the kinetic law distribution vs each type of mass transfer.
  
    Args: 
        fileName: str-file name with which the pdf file save to.

    """  
    (_, df_gen_stat, _, df_gen_stat_PR, _, \
    df_table_PR, df_table_PR_per_model) = self.tuple

    #generate the PR two tables
    try:
      df_table_PR_plot = df_table_PR.div(df_table_PR.sum().sum())
      df_table_PR_per_model_plot = df_table_PR_per_model.div(df_table_PR_per_model.sum().sum())
    except Exception as e:
      raise Exception(e)

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
    fig = plt.figure(figsize = (16,16))
    axes = fig.subplots(nrows=4, ncols=4)
    for i in range(16):
      df_gen_stat_PR_plot[i] = pd.DataFrame(columns = df_gen_stat_PR.columns.tolist())
      df_temp = df_gen_stat_PR[types*i:types*(i+1)]   
      df_gen_stat_PR_plot[i] = pd.concat([df_gen_stat_PR_plot[i],df_temp], ignore_index=True)

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
    fig.savefig(fileName)

  def plotRxnDistOfEachMassTransfer(self, fileName = 'ReactionDistributionOfMassTransfer.pdf'):
    """
    Plot the reaction distribution of each mass transfer.
  
    Args: 
        fileName: str-file name with which the pdf file save to.

    """  
    (_, _, _, _, _, df_table_PR, _) = self.tuple

    #generate the PR two tables
    try:
      df_table_PR_plot = df_table_PR.div(df_table_PR.sum().sum())
    except Exception as e:
      raise Exception(e)

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
    plt.savefig(fileName, dpi=350)

  def plotRxnDistPerModelOfEachMassTransfer(self, fileName = 'ReactionDistributionPerModelOfMassTransfer.pdf'):
    """
    Plot the reaction distribution per model of each mass transfer.
      
    Args: 
        fileName: str-file name with which the pdf file save to.

    """  
    (_, _, _, _, _, _, df_table_PR_per_model) = self.tuple

    #generate the PR two tables
    try:
      df_table_PR_per_model_plot = df_table_PR_per_model.div(df_table_PR_per_model.sum().sum())
    except Exception as e:
      raise Exception(e)

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
    plt.savefig(fileName, dpi=350)

  def _printBriefStatOfKineticLawDistribution(self):
    """
    Print the brief statistics for the kinetic law distribution.
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

  def _printReactionKineticsTypes(self):
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
  #print(analyzer.getKineticLawDistribution(path = 'D:/summer-2020/Jo/', fileName = "KineticLawDistribution.xlsx")) 
  print(analyzer.getKineticLawDistribution(fileName = "KineticLawDistribution.xlsx")) 
  analyzer.plotKineticLawDistribution() 
  print(analyzer.TopFrequentKineticLawType())
  

  print(analyzer.getKineticLawDistributionPerMassTransfer(rct_num=1,prd_num=1,
  fileName="KineticLawDistributionPerMassTransfer.xlsx"))
  analyzer.plotKineticLawDistributionPerMassTransfer(rct_num=1,prd_num=1)
  print(analyzer.TopFrequentKineticLawTypePerMassTransfer(rct_num=1,prd_num=1))
  analyzer.plotKineticLawDistributionVsMassTransfer()

  analyzer.plotRxnDistOfEachMassTransfer()
  analyzer.plotRxnDistPerModelOfEachMassTransfer()

  analyzer._printBriefStatOfKineticLawDistribution()
  analyzer._printReactionKineticsTypes()
  analyzer._saveAllStatisticsInfoToExcel(fileName='statistics_result.xlsx')
 

  print("--- %s seconds ---" % (time.time() - start_time))
