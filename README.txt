============================================

  PASSIVE: Power lAw kineticS iS InjectiVE

============================================

MATLAB was used to develop the function used here.
The function PLKinjective returns whether the chemical reaction network (CRN) is injective or not. Furthermore, the output variables 'injective', 'det_Mstar', and 'model' allow the user to view the following, respectively:

   - 1 or 0 if injective or not, respectively
   - Determinant of the system
   - Complete network with all the species listed in the 'species' field of the structure 'model'

An alternative to PLKinjective is PLKinjectiveLU. This function is used whenever PLKinjective finds it difficult to compute the determinant or to analyze the signs of the determinanat. PLKinjectiveLU computes the diagonal entries of the upper triangle matrix from the LU-decomposition of a chemical reaction network (CRN) for the purpose of computing the determinant to determine whether the network is injective or not. The output variables 'u' and 'model' allow the user to view the following, respectively:

   - Diagonal entries of the upper triangle matrix from the LU-decomposition
   - Complete network with all the species listed in the 'species' field of the structure 'model'



====
Note
====

The idea for the first part of the algorithm comes from [3].



=================================
How to fill out 'model' structure
=================================

'model' is the input for the function conservationLaw. It is a structure, representing the CRN, with the following fields:

   - id: name of the model
   - species: a list of all species in the network; this is left blank since incorporated into the function is a step which compiles all species used in the model
   - reaction: a list of all reactions in the network, each with the following subfields:
        - id: a string representing the reaction
        - reactant: has the following further subfields:
             - species: a list of strings representing the species in the reactant complex
             - stoichiometry: a list of numbers representing the stoichiometric coefficient of each species in the reactant complex (listed in the same order of the species)
        - product: has the following further subfields:
             - species: a list of strings representing the species in the product complex
             - stoichiometry: a list of numbers representing the stoichiometric coefficient of each species in the product complex (listed in the same order of the species)
        - reversible: has the value true or false indicating if the reaction is reversible or not, respectively
        - kinetic: has the following further subfields:
             - reactant1: a list of numbers representing the kinetic order of each species in the reactant complex in the left to right direction (listed in the same order of the species)
             - reactant2: a list of numbers representing the kinetic order of each species in the reactant complex in the right to left direction (listed in the same order of the species) (empty if the reaction is not reversible)

To fill out the 'model' structure, write a string for 'model.id': this is just to put a name to the network. To add the reactions to the network, use the function addReaction where the output is 'model'. addReaction is developed to make the input of reactions of the CRN easier than the input in [7]:

   addReaction
      - OUTPUT: Returns a structure called 'model' with added field 'reaction' with subfields 'id', 'reactant', 'product', 'reversible', and 'kinetic'. The output variable 'model' allows the user to view the network with the added reaction.
      - INPUTS
           - model: a structure, representing the CRN
           - id: visual representation of the reaction, e.g., reactant -> product (string)
           - reactant_species: species of the reactant complex (cell)
           - reactant_stoichiometry: stoichiometry of the species of the reactant complex (cell)
           - reactant_kinetic: kinetic orders of the species of the reactant complex (array)
           - product_species: species of the product complex (cell)
           - product_stoichiometry: stoichiometry of the species of the product complex (cell)
           - product_kinetic: "kinetic orders" of the species of the product complex, if the reaction is reversible (array); if the reaction in NOT reversible, leave blank
           - reversible: logical; whether the reaction is reversible or not (true or false)
      * Make sure the function addReaction is in the same folder/path being used as the current working directory.



========
Examples
========

2 examples are included in this folder:

   - Example 1: Equation 2 in [1]

   - Example 2: Example 1 in [2]



===================
Contact Information
===================

For questions, comments, and suggestions, feel free to contact me at pvnlubenia@yahoo.co.uk.


- Patrick Lubenia (19 March 2024)



==========
References
==========

   [1] Feliu E, Wiuf C (2013) A computational method to preclude multistationarity in networks of interacting species. Bioinformatics 29(18):2327--2334. https://doi.org/10.1093/bioinformatics/btt400

   [2] Hernandez B, Mendoza E, de los Reyes V A (2020) A computational approach to multistationarity of power-law kinetics. J Math Chem 58:56-87. https://doi.org/10.1007/s10910-019-01072-7

   [3] Soranzo N, Altafini C (2009) ERNEST: a toolbox for chemical reaction network theory. Bioinform 25(21):2853-2854. https://doi.org/10.1093/bioinformatics/btp513