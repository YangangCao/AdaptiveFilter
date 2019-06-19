function [ladderVector,...
          kappaVector,...
          posterioriErrorMatrix] = LRLS_pos_ErrorFeedback(desired,input,S)

%   LRLS_pos_ErrorFeedback.m
%       Implements the Error-Feedback LRLS algorithm based on a posteriori errors.
%       (Algorithm 7.3 - book: Adaptive Filtering: Algorithms and Practical
%                                                       Implementation, Diniz)
% 
%   Syntax:
%       [ladderVector,kappaVector,posterioriErrorMatrix] = LRLS_pos_ErrorFeedback(desired,input,S)
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
xiMin_f_Curr             = zeros(1, S.nSectionsLattice +2);
xiMin_f_Prev             = zeros(1, S.nSectionsLattice +2);
xiMin_b_Curr             = zeros(1, S.nSectionsLattice +2);
xiMin_b_Prev             = zeros(1, S.nSectionsLattice +2);
xiMin_b_PPrev            = zeros(1, S.nSectionsLattice +2);
error_b_Curr             = zeros(1, S.nSectionsLattice +2);
error_b_Prev             = zeros(1, S.nSectionsLattice +2);
error_f                  = zeros(1, S.nSectionsLattice +2);
gammaVector_Curr         = zeros(1, S.nSectionsLattice +2);
gammaVector_Prev         = zeros(1, S.nSectionsLattice +2);
kappa_f                  = zeros(1, S.nSectionsLattice +1);
kappa_b                  = zeros(1, S.nSectionsLattice +1);
ladderVector             = zeros(nCoefficients        , nIterations); 
kappaVector              = zeros(nCoefficients        , nIterations); 
posterioriErrorMatrix    = zeros(S.nSectionsLattice +2, nIterations);


%################################################
% Computations
%################################################

%%Initialize Parameters
%
kappa_f                  = zeros(1, S.nSectionsLattice +1);
kappa_b                  = zeros(1, S.nSectionsLattice +1);
deltaVector              = zeros(1, S.nSectionsLattice +1);
%
gammaVector_Prev         = ones(1,S.nSectionsLattice +2);
%
xiMin_f_Prev             = repmat(S.epsilon,1, S.nSectionsLattice +2);
xiMin_b_Prev             = repmat(S.epsilon,1, S.nSectionsLattice +2);
xiMin_b_PPrev            = repmat(S.epsilon,1, S.nSectionsLattice +2); 
%
error_b_Prev             = zeros(1,S.nSectionsLattice +2); 


%-------------------------------------------------------------
%||||||||||||||||||||||| - Main Loop - |||||||||||||||||||||||
%-------------------------------------------------------------

for it = 1:nIterations

    % Set Values for Section 0(Zero)
    gammaVector_Curr(1)         = 1;
    error_b_Curr(1)             = input(it);
    error_f(1)                  = input(it);
    xiMin_f_Curr(1)             = input(it)^2 + S.lambda*xiMin_f_Prev(1);
    xiMin_b_Curr(1)             = xiMin_f_Curr(1);
    posterioriErrorMatrix(1,it) = desired(it);
    
    % Propagate the Order Update Equations
    for ot = 1:S.nSectionsLattice+1
    
        %Delta Time Update
        deltaVector(ot) = S.lambda*deltaVector(ot) + ...
                          (error_b_Prev(ot)*error_f(ot))/gammaVector_Prev(ot);
                      
        %Order Update Equations                        
        gammaVector_Curr(ot+1) = gammaVector_Curr(ot) - ...
                                 (error_b_Curr(ot)^2)/xiMin_b_Curr(ot);
        
        kappa_f(ot) = (gammaVector_Prev(ot+1)/gammaVector_Prev(ot))*...
	              (kappa_f(ot)  + ...
		        (error_b_Prev(ot)*error_f(ot))/...
		        (gammaVector_Prev(ot)*S.lambda*xiMin_b_PPrev(ot)));

        kappa_b(ot) = (gammaVector_Curr(ot+1)/gammaVector_Prev(ot))*...
	              (kappa_b(ot)  + ...
		        (error_b_Prev(ot)*error_f(ot))/...
		        (gammaVector_Prev(ot)*S.lambda*xiMin_f_Prev(ot)));

        error_b_Curr(ot+1) = error_b_Prev(ot) -kappa_b(ot)*error_f(ot);
        error_f(ot+1)      = error_f(ot)      -kappa_f(ot)*error_b_Prev(ot);
        
        xiMin_f_Curr(ot+1) = xiMin_f_Curr(ot) -(deltaVector(ot)^2)/xiMin_b_Prev(ot);
	    xiMin_b_Curr(ot+1) = xiMin_b_Prev(ot) -(deltaVector(ot)^2)/xiMin_f_Curr(ot);
        
        % Feedforward Filtering    
                
        ladderVector(ot,it+1) = (gammaVector_Curr(ot+1)/gammaVector_Curr(ot))*...
	                             (ladderVector(ot,it) + ...
				        (posterioriErrorMatrix(ot,it)*error_b_Curr(ot))/...   
                                        (gammaVector_Curr(ot)*S.lambda*xiMin_b_Prev(ot)));
					
        posterioriErrorMatrix(ot+1,it)=posterioriErrorMatrix(ot,it) - ...
                                      ladderVector(ot,it+1)*error_b_Curr(ot);

        kappaVector(ot,it)      = kappa_f(ot);
 
    end

    
    %Upadate K-1 Info
    gammaVector_Prev = gammaVector_Curr;
    error_b_Prev     = error_b_Curr;
    xiMin_b_PPrev    = xiMin_b_Prev;
    xiMin_b_Prev     = xiMin_b_Curr;
    xiMin_f_Prev     = xiMin_f_Curr;
      
    


end % END OF LOOP

%EOF

