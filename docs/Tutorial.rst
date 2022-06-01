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

   #Statistics 
   print(analyzer.getKTypeDistribution()) 
   print(analyzer.getKTypeDistributionPerMType(rct_num=1,prd_num=1))
   print(analyzer.getMTypeDistribution())
   print(analyzer.getMTypeDistributionPerModel())
   print(analyzer.getNumBiomodelsAnalyzed())
   print(analyzer.getNumRxnsAnalyzed())

   #Query 
   print(analyzer.getTopKType())
   print(analyzer.getKTypeProb(K_type = "NA"))
   print(analyzer.getTopKTypePerMType(rct_num=1,prd_num=1))
   print(analyzer.getKTypeProbPerMType(rct_num=1, prd_num = 1, K_type="NA"))
   print(analyzer.getTopMType())
   print(analyzer.getMTypeProb(rct_num = 1, prd_num = 1))

   #Presentations
   #plot
   analyzer.plotKTypeDistribution(path = 'D:/path/to/folder/')
   analyzer.plotKTypeDistributionPerMType(rct_num=1,prd_num=1)
   analyzer.plotKTypeDistributionVsMType()
   analyzer.plotMtypeDistribution()
   analyzer.plotMTypeDistributionPerModel()
   #table
   analyzer.tableKTypeDistribution(path = 'D:/path/to/folder/', fileName = "Kinetics.xlsx")
   analyzer.tableKTypeDistribution()
   analyzer.tableKTypeDistributionPerMType(rct_num=1,prd_num=1)
   analyzer.tableMTypeDistribution()
   analyzer.tableMTypeDistributionPerModel()