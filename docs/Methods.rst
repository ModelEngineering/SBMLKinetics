.. _Methods:
 

Methods
=============

--------------------------------
Define kinetics type (K type)
--------------------------------

- Zeroth order: the number of species in the kinetic law is zero;
- Uni-directional mass action: the kinetic law is a single product of terms and all species in the kinetic law are reactants;
- Uni-term with the moderator: the kinetic law is a single product of terms and NOT UNDR;
- Bi-directional mass action: the kinetic law is the difference between two products of terms and the first product of terms with all reactants while the second product of terms with all products;
- Bi-terms with the moderator: the kinetic law is the difference between two products of terms and NOT BIDR;
- Michaelis-Menten kinetics without explicit enzyme: the kinetic law is in the format of Michaelis-Menten expressions without an explicit enzyme;
- Michaelis-Menten kinetics with an explicit enzyme: the kinetic law is in the format of Michaelis-Menten expressions with an explicit enzyme;
- Hill equation: the kinetic law is in the format of Hill equations;
- Fraction format other than MM, MMCAT and HILL: the kinetic law is in the format of fraction with at least one species in the denominator and NOT MM, MMCAT or HILL;
- Not classified: not classified kinetics.

To classify the K types, we have considered the following features as kinetics properties (K properties).

- a. the number of species in the kinetic law;
- b. whether the kinetic law is a single product of terms;
- c. whether the kinetic law is the difference between two products of terms;
- d. whether the first (/second) product of terms with all reactants (/products);
- e. whether the kinetic law is in the format of Michaelis-Menten expressions with the single reactant in the numerator;
- f. whether the kinetic law is in the format of Michaelis-Menten expressions with the product of the single reactant and another species in the numerator;
- g. whether the kinetic law is in the format of Hill equation;
- h. whether the kinetic law is in the fraction format with at least one species in the denominator.

For example, ZERO means there are no species in the kinetic law. Therefore, the number of species 
in the kinetic law should be zero referring to a. In addition, the kinetics of UNDR is a 
single product of terms referring to b and the species in the product of terms are all 
reactants referring to d. See Table below for details about how K properties determine the K types
in practice. 

.. list-table:: 
   :widths: 20 20 20 20 20 20 20 20 20 20
   :header-rows: 1

   * - K properties
     - ZERO
     - UNDR
     - UNMO
     - BIDR
     - BIMO
     - MM
     - MMCAT
     - HILL
     - FR
   * - a.
     - = 0
     - > 0
     - > 0
     - > 0
     - > 0
     - 1
     - 2
     - 1
     - > 0
   * - b.
     - 
     - Yes
     - Yes
     - 
     - 
     - 
     - 
     - 
     -
   * - c.
     - 
     - 
     - 
     - Yes
     - Yes
     - 
     - 
     - 
     -
   * - d.
     - 
     - Yes
     - No
     - Yes
     - No
     - 
     - 
     - 
     -
   * - e.
     - 
     - 
     - 
     - 
     - 
     - Yes
     - 
     - 
     - No
   * - f.
     - 
     - 
     - 
     - 
     - 
     - 
     - Yes
     - 
     - No
   * - g.
     - 
     - 
     - 
     - 
     - 
     - 
     - 
     - Yes
     - No
   * - h.
     - 
     - 
     - 
     - 
     - 
     - Yes
     - Yes
     - Yes
     - Yes


--------------------------------
Define reaction type (R type)
--------------------------------

R type is defined mainly by the structure of reactions. It is a simplification of reaction 
stoichiometry considering the number of reactants and the number of products. For example, 
there are two reactants and one product in the reaction of A + B -> C. 
We define R type quantitatively with a pair of the number of reactants (R = 0, 1, 2, >2) and 
products (P = 0, 1, 2, >2). As a comparison, the stoichiometry matrix indicates reactants and 
products by negative and positive signs respectively.

