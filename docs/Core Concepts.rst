.. _Core Concepts:
 

Core Concepts
=============
-------
K Type
-------
Kinetics type (K type) including ten types:

- "ZERO" (Zeroth order);
- "UNDR" (Uni-directional mass action);
- "UNMO" (Uni-term with the moderator);
- "BIDR" (Bi-directional mass action);
- "BIMO" (Bi-terms with the moderator);
- "MM" (Michaelis-Menten kinetics without an explicit enzyme);
- "MMCAT" (Michaelis-Menten kinetics with an explicit enzyme);
- "HILL" (Hill equations);
- "FR" (Kinetics in the format of fraction other than MM, MMCAT or HILL);
- "NA" (not classified kinetics). 

Example: K_type = SBMLKinetics.types.K_type("NA").


-------
R Type
-------
Reaction type (R type) is quantitatively represented by the number of reactants 
(R = 0, 1, 2, 3 (representing>2)) and products (P = 0, 1, 2, 3 (representing>2)).

Example: R_Type = SBMLKinetics.types.R_Type(1,1).

