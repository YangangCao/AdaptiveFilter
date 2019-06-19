function [ladderVector,...
          kappaVector,...
          posterioriErrorMatrix] = NLRLS_pos(desired,input,S)

%   NLRLS_pos.m
%       Implements the Normalized Lattice RLS algorithm based on a posteriori error.
%       (Algorithm 7.2 - book: Adaptive Filtering: Algorithms and Practical
%                                                       Implementation, Diniz)
% 
%   Syntax:
%       [ladderVector,kappaVector,posterioriErrorMatrix] = NLRLS_pos(desired,input,S)
% 
%   Input Arguments: 
%       . desired           : Desired Signal.                          (ROW Vector)
%       . input             : Signal fed into the adaptive filter.     (ROW Vector)
%       . S                 : Structure with the following fields
%           - lambda             : Forgetting factor.                  (0 << lambda < 1)
%           - nSectionsLattice   : Number of lattice sections, refered as N in the textbook.
%           - epsilon            : Initilization of xiMin_backward and xiMin_forward. 
% 
%   Output Arguments:
%       . ladderVector          : Store the ladder coefficients for each iteration, refered 
%                                 as v in the textbook.
%                                 (Coefficients at one iteration are column vector)
%       . kappaVector           : Store the reflection coefficients for each iteration, considering 
%                                 only the forward part.
%                                 (Coefficients at one iteration are column vector)
%       . posterioriErrorMatrix : Store the a posteriori errors at each section of the lattice for 
%                                 each iteration.
%                                 (The errors at one iteration are column vectors)
%
%   Authors:
%       . Guilherme de Oliveira Pinto   - guilhermepinto7@gmail.com  &  guilherme@lps.ufrj.br
%       . Markus Vinícius Santos Lima   - mvsl20@gmailcom            &  markus@lps.ufrj.br
%       . Wallace Alves Martins         - wallace.wam@gmail.com      &  wallace@lps.ufrj.br
%       . Luiz Wagner Pereira Biscainho - cpneqs@gmail.com           &  wagner@lps.ufrj.br
%       . Paulo Sergio Ramirez Diniz    -                               diniz@lps.ufrj.br
%



%################################################
% Data Initialization
%################################################

% Basic Procedure
nCoefficients           = S.nSectionsLattice +1;  
nIterations             = length(desired); 

% Pre Allocations
deltaVector              = zeros(1, S.nSectionsLattice +1);
deltaVector_D            = zeros(1, S.nSectionsLattice +1);   %acho que deveria ser +2
error_b_Curr             = zeros(1, S.nSectionsLattice +2);
error_b_Prev             = zeros(1, S.nSectionsLattice +2);
error_f                  = zeros(1, S.nSectionsLattice +2);
ladderVector             = zeros(nCoefficients        , nIterations); 
kappaVector              = zeros(nCoefficients        , nIterations); 
posterioriErrorMatrix    = zeros(S.nSectionsLattice +2, nIterations);


%################################################
% Computations
%################################################

%%Initialize Parameters
%
deltaVector              = zeros(1, S.nSectionsLattice+1);
deltaVector_D            = zeros(1, S.nSectionsLattice+1);
%
error_b_Prev             = zeros(1,S.nSectionsLattice +2); 
%
sigma_d_2                = S.epsilon;
sigma_x_2                = S.epsilon;


%-------------------------------------------------------------
%||||||||||||||||||||||| - Main Loop - |||||||||||||||||||||||
%-------------------------------------------------------------

for it = 1:nIterations

    sigma_x_2                   = S.lambda*sigma_x_2 + input(it)^2;
    sigma_d_2                   = S.lambda*sigma_d_2 + desired(it)^2;
    
    % Set Values for Section 0(Zero)
    error_b_Curr(1)             = input(it)/(sigma_x_2^0.5);
    error_f(1)                  = error_b_Curr(1);
    posterioriErrorMatrix(1,it) = desired(it)/(sigma_d_2^0.5);
    
    % Propagate the Order Update Equations
    for ot = 1:S.nSectionsLattice+1
    
        %Delta Time Update
        deltaVector(ot) = deltaVector(ot)*...
                          sqrt((1 -error_b_Prev(ot)^2)*(1- error_f(ot)^2)) + ...
                          error_b_Prev(ot)*error_f(ot);

        %Order Update Equations                        
        error_b_Curr(ot+1) = (error_b_Prev(ot) -deltaVector(ot)*error_f(ot))/...
	                          sqrt((1 -deltaVector(ot)^2)*(1 -error_f(ot)^2));
        error_f(ot+1)      = (error_f(ot) - deltaVector(ot)*error_b_Prev(ot))/...
                           sqrt((1 -deltaVector(ot)^2)*(1 -error_b_Prev(ot)^2));
 
        % Feedforward Filtering    
        deltaVector_D(ot) = deltaVector_D(ot)*...
                            sqrt((1 -error_b_Curr(ot)^2)*(1 -posterioriErrorMatrix(ot,it)^2)) + ...
                            posterioriErrorMatrix(ot,it)*error_b_Curr(ot);
        
        
posterioriErrorMatrix(ot+1,it)= 1/sqrt((1 -error_b_Curr(ot)^2)*(1 -deltaVector_D(ot)^2))*...
                                (posterioriErrorMatrix(ot,it) -deltaVector_D(ot)*error_b_Curr(ot));

         % Outputs:
         ladderVector(ot,it) = deltaVector_D(ot);
         kappaVector(ot,it)  = deltaVector(ot);
 
    end
    
    %Upadate K-1 Info
    error_b_Prev     = error_b_Curr;


end % END OF LOOP

%EOF

