.. _Tutorial:
 

Tutorial
=============

.. code-block:: python

   import SBMLKinetics

   initial_model_indx = 5
   final_model_indx = 6

   model_indices = range(initial_model_indx, final_model_indx+1)
   analyzer = KineticAnalyzer(path = 'D:/path/to/folder',
   dataSet = "biomodels.zip", model_indices=model_indices) 

   #Query Distributions 
   print(analyzer.getKTypeDistribution()) 
   print(analyzer.getKTypeDistributionPerMType(M_type = SBMLKinetics.types.M_type(1,1)))
   print(analyzer.getMTypeDistribution())
   print(analyzer.getMTypeDistributionPerModel())

   #Query Elements
   print(analyzer.getTopKType()[0].K_type_str)
   print(analyzer.getKTypeProb(K_type = SBMLKinetics.types.K_type("NA")))
   print(analyzer.getTopKTypePerMType(M_type = SBMLKinetics.types.M_type(1,1))[0].K_type_str)
   print(analyzer.getKTypeProbPerMType(M_type = SBMLKinetics.types.M_type(1,1), K_type=SBMLKinetics.types.K_type("NA")))
   print(analyzer.getTopMType()[0].rct_num)
   print(analyzer.getMTypeProb(M_type = SBMLKinetics.types.M_type(1,1)))
   print(analyzer.getNumBiomodelsAnalyzed())
   print(analyzer.getNumRxnsAnalyzed())

   #Plots
   analyzer.plotKTypeDistribution(path = 'D:/path/to/folder/')
   analyzer.plotKTypeDistributionPerMType(M_type = SBMLKinetics.types.M_type(1,1))
   analyzer.plotKTypeDistributionVsMType()
   analyzer.plotMtypeDistribution()
   analyzer.plotMTypeDistributionPerModel()
