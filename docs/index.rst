.. SBMLKinetics documentation master file, created by
   sphinx-quickstart on Mon Nov  8 14:18:50 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to SBMLKinetics's documentation!
========================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   Introduction
   Core Concepts
   Tutorial
   Examples
   Query Methods
   Plot Methods
    
    
------------
Introduction
------------

SBMLKinetics is a Python package to evaluate and classify kinetics expressions in SBML models. 
There are many possible kinetics laws like zeroth order, mass action, Michaelis-Menten, 
Hill equations and others. This work characterizes the kinetic laws used in the BioModels 
Database to improve modeling best practices. Our tool can analyze any dataset with SBML files 
as input. Users can also use this tool to compare different data sets. For instance, we 
compare the distribution of kinetic laws for signaling and metabolic networks and find the 
substantial differences between two types of networks. If you are using any of the code, 
please cite the PYPI web page (https://pypi.org/project/SBMLKinetics/).
