"""
This script is to look for uninitialized variables.
Make sure that you have setup your PYTHONPATH environment
variable as described in the github repository.
"""

# Import the required files
from src.common.simple_sbml import SimpleSBML
import src.common.simple_sbml as simple_sbml
import src.common.constants as cn

import numpy as np
import os

#self add:
from libsbml import *

iterator = simple_sbml.modelIterator(initial=0, final=20)
for idx,item in enumerate(iterator):
  if item is not None:
    name = item.filename
    print(name)

    # Create an SBML model. We'll use the model
    # data/
    
    path = os.path.join(cn.PROJECT_DIR, "data")
    path = os.path.join(path, name)

    simple = item.model # Creates a model
    model = simple.model

    species_num = model.getNumSpecies()
    #print(species_num)
    parameter_num = model.getNumParameters()
    #print(parameter_num)

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

    #print("species:", species_list)
    #print("parameter:", parameter_list)

    uninitialized_species = []
    uninitialized_parameters = []
    for reaction in simple.reactions:
      reactant_list = []
      product_list = []
      species_in_kinetic_law = []
      parameters_in_kinetic_law = []
      reactant_stg = " + ".join(
        [r.getSpecies() for r in reaction.reactants])
      reactant_list.append([r.getSpecies() for r in reaction.reactants])
      product_stg = " + ".join(
        [p.getSpecies() for p in reaction.products])
      product_list.append([p.getSpecies() for p in reaction.products])

      species_parameter_list = list(dict.fromkeys(reaction.kinetic_law.symbols))
      if(len(reactant_list[0]) != 0):
        species_parameter_list.append(reactant_list[0][0])
      if(len(product_list[0]) != 0):
        species_parameter_list.append(product_list[0][0])
      species_parameter_list = list(dict.fromkeys(species_parameter_list))

      for i in range(len(species_parameter_list)):
        if species_parameter_list[i] in species_list:
          species_in_kinetic_law.append(species_parameter_list[i])
        elif species_parameter_list[i] in parameter_list:
          parameters_in_kinetic_law.append(species_parameter_list[i])

     
      for i in range(len(species_in_kinetic_law)):
        species = model.getSpecies(species_in_kinetic_law[i])
        species_id = species.getId()
        species_concentration = species.getInitialConcentration()
        species_bool_concentration = species.isSetInitialConcentration()
        if species_bool_concentration == False:
          uninitialized_species.append(species_id)

      for i in range(len(parameters_in_kinetic_law)):
        parameter = model.getParameter(parameters_in_kinetic_law[i])

        parameter_id =  parameter.getId()
        parameter_bool = parameter.isSetValue()
        if parameter_bool == False:
          uninitialized_parameters.append(parameter_id)

    #remove the repeated recorded uninitiazlied variables (parameters, speci>
    uninitialized_species = list(dict.fromkeys(uninitialized_species))
    uninitialized_parameters = list(dict.fromkeys(uninitialized_parameters))

    uninitialized_species_final = []
    uninitialized_parameters_final = []
    for i in range(len(uninitialized_species)):
      if model.getAssignmentRule(uninitialized_species[i]) == None:
        uninitialized_species_final.append(uninitialized_species[i])
  
    for i in range(len(uninitialized_parameters)):
      if model.getAssignmentRule(uninitialized_parameters[i]) == None:
        uninitialized_parameters_final.append(uninitialized_parameters[i])

    print("uninitialized species: %s" % uninitialized_species_final)
    print("uninitialized parameters: %s" % uninitialized_parameters_final)

