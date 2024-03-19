% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                                           %
%    PLKinjectiveLU                                                         %
%                                                                           %
%                                                                           %
% OUTPUT: Computes the diagonal entries of the upper triangle matrix from   %
%    the LU-decomposition of a chemical reaction network (CRN) for the      %
%    purpose of computing the determinant to determine whether the network  %
%    is injective or not. The output variables 'u' and 'model' allow the    %
%    user to view the following, respectively:                              %
%       - Diagonal entries of the upper triangle matrix from the LU-        %
%            decomposition                                                  %
%       - Complete network with all the species listed in the 'species'     %
%            field of the structure 'model'                                 %
%                                                                           %
% INPUT: model: a structure, representing the CRN (see README.txt for       %
%    details on how to fill out the structure)                              %
%                                                                           %
% Notes:                                                                    %
%    1. It is assumed that the CRN has power law kinetics.                  %
%    2. This version is for big systems who determinant is difficult to     %
%          compute or analysis for signs: hence, u is an output so the user %
%          can opt to manually compute for the determinant (simply the      %
%          product of the entries of u) or investigate its sign.            %
%    3. Ideas for some parts of the code was motivated by [2].              %
%                                                                           %
% References                                                                %
%    [1] ﻿Feliu E, Wiuf C (2013) A computational method to preclude          %
%           multistationarity in networks of interacting species.           %
%           Bioinformatics 29(18):2327--2334.                               %
%           https://doi.org/10.1093/bioinformatics/btt400                   %
%    [2] Soranzo N, Altafini C (2009) ERNEST: a toolbox for chemical        %
%           reaction network theory. Bioinform 25(21):2853–2854.            %
%           https://doi.org/10.1093/bioinformatics/btp513                   %
%                                                                           %
% Created: 21 January 2023                                                  %
% Last Modified: 19 March 2024                                              %
%                                                                           %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



function [u, model] = PLKinjectiveLU(model)

%
% Step 1: Create a list of all species indicated in the reactions
%

% Initialize list of species
model.species = { };

% Get all species from reactants
for i = 1:numel(model.reaction)
    for j = 1:numel(model.reaction(i).reactant)
        model.species{end+1} = model.reaction(i).reactant(j).species;
    end
end

% Get species from products
for i = 1:numel(model.reaction)
    for j = 1:numel(model.reaction(i).product)
        model.species{end+1} = model.reaction(i).product(j).species;
    end
end

% Get only unique species
model.species = unique(model.species);

% Count the number of species
m = numel(model.species);



%
% Step 2: Form stoichiometric matrix N
%

% Initialize the matrix of reactant complexes
reactant_complex = [ ];

% Initialize the matrix of product complexes
product_complex = [ ];

% Initialize the stoichiometric matrix
N = [ ];

% For each reaction in the model
for i = 1:numel(model.reaction)
  
    % Initialize the vector for the reaction's reactant complex
    reactant_complex(:, end+1) = zeros(m, 1);
    
    % Fill it out with the stoichiometric coefficients of the species in the reactant complex
    for j = 1:numel(model.reaction(i).reactant)
        reactant_complex(find(strcmp(model.reaction(i).reactant(j).species, model.species), 1), end) = model.reaction(i).reactant(j).stoichiometry;
    end
    
    % Initialize the vector for the reaction's product complex
    product_complex(:, end+1) = zeros(m, 1);
    
    % Fill it out with the stoichiometric coefficients of the species in the product complex
    for j = 1:numel(model.reaction(i).product)
        product_complex(find(strcmp(model.reaction(i).product(j).species, model.species), 1), end) = model.reaction(i).product(j).stoichiometry;
    end
    
    % Create a vector for the stoichiometric matrix: Difference between the two previous vectors
    N(:, end+1) = product_complex(:, end) - reactant_complex(:, end);
    
    % If the reaction is reversible
    if model.reaction(i).reversible
      
        % Insert a new vector for the reactant complex: make it same as the product complex
        reactant_complex(:, end+1) = product_complex(:, end);
        
        % Insert a new vector for the product complex: make it the same as the reactant complex
        product_complex(:, end+1) = reactant_complex(:, end-1);
        
        % Insert a new vector in the stoichiometric matrix: make it the additive inverse of the vector formed earlier
        N(:, end+1) = -N(:, end);
    end
end

% Count the total number of reactions
r = size(N, 2);



%
% Step 3: Form kinetic order matrix F
%    Reminder: Algorithm assumes the system is power law

% Initialize matrix F
F = [ ];

% Go through each reaction
for i = 1:numel(model.reaction)
    
    % Case 1: The reaction is NOT reversible
    if model.reaction(i).reversible == 0
        
        % Add a row of zeros
        F(end+1, :) = zeros(1, m);
        
        % Fill out the kinetic order of all the species in the reactant
        for j = 1:numel(model.reaction(i).kinetic.reactant1)
            F(end, find(strcmp(model.reaction(i).reactant(j).species, model.species))) = model.reaction(i).kinetic.reactant1(j);
        end
    
    % Case 2: The reaction is reversible
    else
        
        % Add a row of zeros
        F(end+1, :) = zeros(1, m);
        
        % Fill out the kinetic order of all the species in the reactant in the first direction
        for j = 1:numel(model.reaction(i).kinetic.reactant1)
            F(end, find(strcmp(model.reaction(i).reactant(j).species, model.species))) = model.reaction(i).kinetic.reactant1(j);
        end
        
        % Add a row of zeros
        F(end+1, :) = zeros(1, m);
        
        % Fill out the kinetic order of all the species in the reactant in the other direction
        for j = 1:numel(model.reaction(i).kinetic.reactant2)
            F(end, find(strcmp(model.reaction(i).product(j).species, model.species))) = model.reaction(i).kinetic.reactant2(j);
        end
    end
end



%
% Step 4: Compute a basis of conservation laws
%

b = rref(null(N', 'r')');



%
% Step 5: Form matrix M
%

% Get the number of complexes
n = size(N, 2);

% Form diagonal matrix with positive entries k
K = diag(sym(strcat('k', string(1:m))));
assume(diag(K) > 0);

% Form diagonal matrix with positive entries z
Z = diag(sym(strcat('z', string(1:r))));
assume(diag(Z) > 0);

% Form matrix M
% N: mxr
% Z: rxr
% F: rxm
% K: mxm
M = N*Z*F*K;



%
% Step 6: Form matrix Mstar
%

% Initialize Mstar
Mstar = M;

% Determine the rank of the CRN
s = rank(N);

% Get the index of the first nonzero element of each basis
[~, nonzero_index] = max(b ~= 0, [ ], 2);

% Replace the row of M based on the first nonzero element of each basis
if s < m
    for i = 1:size(nonzero_index)
        Mstar(nonzero_index(i),:) = b(i, :);
    end
end



%
% Step 7: Compute the determinant of Mstar
%

% det_Mstar = det(Mstar)

% ALTERNATIVE METHOD

% Perform LU factorization on Mstar
[L, U, P] = lu(Mstar);

% Get the diagonal entries of U
u = diag(U);

% Output message
fprintf('\nYou may now manually compute for the determinant by taking the product of the entries of u, then determine the signs of this determinant.\n\n')

% % Compute the determinant using the upper triangular matrix
% det_Mstar = expand(prod(u));
% 
% 
% 
% %
% % Step 8: Check if The polynomial determinant has all positive or all negative terms
% %
% 
% % Get the coefficients of the polynomial determinant
% coefficients = coeffs(det_Mstar);
% 
% % If the determinant is the zero polynomial
% if det_Mstar == 0
%     injective = 0;
%     fprintf('\nThe network is NOT injective since the determinant is identically zero.\n\n')
% else
% 
%     % If all coefficients of the polynomial are positive or all negative
%     if or(all(coefficients > 0), all(coefficients < 0))
%         injective = 1;
%         fprintf('\nThe network is injective.\n\n')
%     else
%         injective = 0;
%         fprintf('\nThe network is NOT injective.\n\n')
%     end
% end

end