function [posterioriErrorVector,...
          prioriErrorVector,...
          coefficientVector] = Fast_RLS(desired,input,S)

%   Fast_RLS.m
%       Implements the Fast Transversal RLS algorithm for REAL valued data.
%       (Algorithm 8.1 - book: Adaptive Filtering: Algorithms and Practical
%                                                       Implementation, 3rd Ed., Diniz)
% 
%   Syntax:
%       [posterioriErrorVector,prioriErrorVector,coefficientVector] = Fast_RLS(desired,input,S)
% 
%   Input Arguments: 
%       . desired           : Desired Signal.                          (ROW Vector)
%       . input             : Signal fed into the adaptive filter.     (ROW Vector)
%       . S                 : Structure with the following fields
%           - lambda             : Forgetting factor.                  (0 << lambda < 1)  
%           - predictorOrder     : refered as N in the textbook.
%           - epsilon            : Initilization of xiMin_backward and xiMin_forward. (usually 0 < epsilon <= 1)
% 
%   Output Arguments:
%       . posterioriErrorVector  : Store the a posteriori error for each iteration.      (COLUMN vector)
%       . prioriErrorVector      : Store the a priori error for each iteration.          (COLUMN vector)
%       . coefficientVector      : Store the estimated coefficients for each iteration.
%                                  (Coefficients at one iteration are COLUMN vector)
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

% Basic Parameters
nCoefficients           = S.predictorOrder +1;  
nIterations             = length(desired); 

% Pre Allocations
xiMin_f_Curr             = 0;
xiMin_f_Prev             = 0;
xiMin_b                  = 0;
gamma_Np1                = 0;
gamma_N                  = 0;
w_f                      = zeros(nCoefficients, 1);
w_b                      = zeros(nCoefficients, 1);
coefficientVector        = zeros(nCoefficients, nIterations+1); 
error_f                  = 0;
error_f_line             = 0;
error_b                  = 0;
error_b_line             = 0;
posterioriErrorVector    = zeros(nIterations,1);
prioriErrorVector        = zeros(nIterations,1);
phiHatN                  = zeros(nCoefficients,1);  
phiHatNPlus1             = zeros(nCoefficients+1,1);
regressor                = zeros(nCoefficients+1,1);


%################################################
% Computations
%################################################

%%Initialize Parameters
w_f                      = zeros(nCoefficients, 1);
w_b                      = zeros(nCoefficients, 1);
coefficientVector(:,1)   = zeros(nCoefficients, 1); 
%
phiHatN                  = zeros(nCoefficients, 1); 
%
gamma_N                  = 1;
%
xiMin_f_Prev             = S.epsilon;
xiMin_b                  = S.epsilon;
%
prefixedInput            = [ zeros(1,nCoefficients) input];


%-------------------------------------------------------------
%||||||||||||||||||||||| - Main Loop - |||||||||||||||||||||||
%------------------------------------------------------------

for it = 1:nIterations

    regressor    = prefixedInput(it+nCoefficients:-1:it).';     
    
    error_f_line = regressor.'*[ 1 ;  -w_f];
    error_f      = error_f_line*gamma_N;
    xiMin_f_Curr = S.lambda*xiMin_f_Prev + error_f_line*error_f;
    phiHatNPlus1 = [0 ; phiHatN ]  + ...
                   1/(S.lambda*xiMin_f_Prev)*[1 ; -w_f]*error_f_line;
    w_f          = w_f + phiHatN*error_f;
    
    gamma_Np1    = (S.lambda*xiMin_f_Prev*gamma_N)/xiMin_f_Curr;
    error_b_line = S.lambda*xiMin_b*phiHatNPlus1(end);
    gamma_N      = 1/(1/gamma_Np1 -phiHatNPlus1(end)*error_b_line);
    
    error_b      = error_b_line*gamma_N;
    xiMin_b      = S.lambda*xiMin_b + error_b*error_b_line;
    phiHatN      = phiHatNPlus1(1:end-1) +phiHatNPlus1(end)*w_b; 
    w_b          = w_b  +  phiHatN*error_b;  
    
    
    % Joint Process Estimation
    prioriErrorVector(it)     = desired(it) -...
                                coefficientVector(:,it).'*regressor(1:end-1);
    posterioriErrorVector(it) = prioriErrorVector(it)*gamma_N;
    coefficientVector(:,it+1) = coefficientVector(:,it) + ...
                                phiHatN*posterioriErrorVector(it);
    
    % K- 1 Info
    xiMin_f_Prev     = xiMin_f_Curr;


end % END OF LOOP


%EOF

