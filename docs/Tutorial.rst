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
   analyzer.getKTypeDistribution()
   analyzer.getKTypeDistributionPerMType(M_type = SBMLKinetics.types.M_type(1,1))
   analyzer.getMTypeDistribution()
   analyzer.getMTypeDistributionPerModel()

   #Query Elements
   analyzer.getTopKType()[0].K_type_str
   analyzer.getKTypeProb(K_type = SBMLKinetics.types.K_type("NA")
   analyzer.getTopKTypePerMType(M_type = SBMLKinetics.types.M_type(1,1))[0].K_type_str
   analyzer.getKTypeProbPerMType(M_type = SBMLKinetics.types.M_type(1,1), K_type=SBMLKinetics.types.K_type("NA"))
   analyzer.getTopMType()[0].rct_num
   analyzer.getMTypeProb(M_type = SBMLKinetics.types.M_type(1,1))
   analyzer.getNumBiomodelsAnalyzed()
   analyzer.getNumRxnsAnalyzed()

   #Plots
   analyzer.plotKTypeDistribution(path = 'D:/path/to/folder/')
   analyzer.plotKTypeDistributionPerMType(M_type = SBMLKinetics.types.M_type(1,1))
   analyzer.plotKTypeDistributionVsMType()
   analyzer.plotMtypeDistribution()
   analyzer.plotMTypeDistributionPerModel()

For example, the output of analyzer.getKTypeDistribution() is a dataframe as below.

.. list-table:: 
   :widths: 25 25 25 25 25
   :header-rows: 1

   * - Classifications
     - Percentage
     - Percentage standard error
     - Percentage per model
     - Percentage per model standard error
   * - ZERO
     - 0.33333
     - 0
     - 0.33333
     - 0
   * - UNDR
     - 0.33333
     - 0
     - 0.33333
     - 0
   * - UNMO
     - 0
     - 0
     - 0
     - 0
   * - BIDR
     - 0
     - 0
     - 0
     - 0
   * - BIMO
     - 0
     - 0
     - 0
     - 0
   * - MM
     - 0
     - 0
     - 0
     - 0
   * - MMCAT
     - 0
     - 0
     - 0
     - 0
   * - HILL
     - 0
     - 0
     - 0
     - 0
   * - FR
     - 0
     - 0
     - 0
     - 0
   * - NA
     - 0.33333
     - 0
     - 0.33333
     - 0
