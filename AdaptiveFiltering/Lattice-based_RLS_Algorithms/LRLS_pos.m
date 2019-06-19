function [ladderVector,...
          kappaVector,...
          posterioriErrorMatrix] = LRLS_pos(desired,input,S)

%   LRLS_pos.m
%       Implements the Lattice RLS algorithm based on a posteriori errors.
%       (Algorithm 7.1 - book: Adaptive Filtering: Algorithms and Practical
%                                                       Implementation, Diniz)
% 
%   Syntax:
%       [ladderVector,kappaVector,posterioriErrorMatrix] = LRLS_pos(desired,input,S)
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
deltaVector_D            = zeros(1, S.nSectionsLattice +1);
xiMin_f                  = zeros(1, S.nSectionsLattice +2);
xiMin_b_Curr             = zeros(1, S.nSectionsLattice +2);
xiMin_b_Prev             = zeros(1, S.nSectionsLattice +2);
gammaVector_Curr         = zeros(1, S.nSectionsLattice +2);
gammaVector_Prev         = zeros(1, S.nSectionsLattice +2);
error_b_Curr             = zeros(1, S.nSectionsLattice +2);
error_b_Prev             = zeros(1, S.nSectionsLattice +2);
error_f                  = zeros(1, S.nSectionsLattice +2);
ladderVector             = zeros(nCoefficients        , nIterations); 
kappaVector              = zeros(nCoefficients        , nIterations); 
posterioriErrorMatrix    = zeros(S.nSectionsLattice +2, nIterations);
%  backwardErrorMatrix      = zeros(S.nSectionsLattice +2, nIterations);
kappa_f                  = 0;
kappa_b                  = 0;


%################################################
% Computations
%################################################

%%Initialize Parameters
%
deltaVector              = zeros(1, S.nSectionsLattice+1);
deltaVector_D            = zeros(1, S.nSectionsLattice+1);
%
xiMin_f                  = repmat(S.epsilon,1, S.nSectionsLattice +2);
xiMin_b_Prev             = repmat(S.epsilon,1, S.nSectionsLattice +2);
%
gammaVector_Prev         =  ones(1,S.nSectionsLattice +2);
error_b_Prev             = zeros(1,S.nSectionsLattice +2); 


%-------------------------------------------------------------
%||||||||||||||||||||||| - Main Loop - |||||||||||||||||||||||
%-------------------------------------------------------------

for it = 1:nIterations

    % Set Values for Section 0(Zero)
    gammaVector_Curr(1)         = 1;
    error_b_Curr(1)             = input(it);
    error_f(1)                  = input(it);
    xiMin_f(1)                  = input(it)^2 + S.lambda*xiMin_f(1);
    xiMin_b_Curr(1)             = xiMin_f(1);
    posterioriErrorMatrix(1,it) = desired(it);
    
    % Propagate the Order Update Equations
    for ot = 1:S.nSectionsLattice+1
    
        %Delta Time Update
        deltaVector(ot) = S.lambda*deltaVector(ot) + ...
                          (error_b_Prev(ot)/gammaVector_Prev(ot))*error_f(ot);
                      
        %Order Update Equations                        
        gammaVector_Curr(ot+1) = gammaVector_Curr(ot) - ...
                            (error_b_Curr(ot)^2)/xiMin_b_Curr(ot);
        
        kappa_b = deltaVector(ot)/xiMin_f(ot);
        kappa_f = deltaVector(ot)/xiMin_b_Prev(ot);
        
        error_b_Curr(ot+1) = error_b_Prev(ot) -kappa_b*error_f(ot);
        error_f(ot+1)      = error_f(ot)      -kappa_f*error_b_Prev(ot);
        
        xiMin_b_Curr(ot+1) = xiMin_b_Prev(ot) -deltaVector(ot)*kappa_b; 
        xiMin_f(ot+1)      = xiMin_f(ot)      -deltaVector(ot)*kappa_f;
        
        % Feedforward Filtering    
        %Delta_D Time Update
        deltaVector_D(ot) = S.lambda*deltaVector_D(ot) + ...
                            (error_b_Curr(ot)/gammaVector_Curr(ot))*...
                             posterioriErrorMatrix(ot,it);
        
        ladderVector(ot,it) = deltaVector_D(ot)/xiMin_b_Curr(ot);
        
        posterioriErrorMatrix(ot+1,it)=posterioriErrorMatrix(ot,it) - ...
                                      ladderVector(ot,it)*error_b_Curr(ot);

        kappaVector(ot,it)      = kappa_f;
 
    end

    
    %backwardErrorMatrix(:,it)  =  error_b_Curr.';
    %Upadate K-1 Info
    gammaVector_Prev = gammaVector_Curr;
    error_b_Prev     = error_b_Curr;
    xiMin_b_Prev     = xiMin_b_Curr;


end % END OF LOOP

%EOF
