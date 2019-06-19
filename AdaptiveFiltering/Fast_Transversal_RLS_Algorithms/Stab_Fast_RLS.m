function [posterioriErrorVector,...
          prioriErrorVector,...
          coefficientVector] = Stab_Fast_RLS(desired,input,S)

%   Stab_Fast_RLS.m
%       Implements the Stabilized Fast Transversal RLS algorithm for REAL valued data.
%       (Algorithm 8.2 - book: Adaptive Filtering: Algorithms and Practical
%                                                       Implementation, 3rd Ed., Diniz)
% 
%   Syntax:
%       [posterioriErrorVector,prioriErrorVector,coefficientVector] = Stab_Fast_RLS(desired,input,S)
% 
%   Input Arguments: 
%       . desired           : Desired Signal.                          (ROW Vector)
%       . input             : Signal fed into the adaptive filter.     (ROW Vector)
%       . S                 : Structure with the following fields
%           - lambda             : Forgetting factor.                  (0 << lambda < 1) 
%           - predictorOrder     : refered as N in the textbook.
%           - epsilon            : Initilization of xiMin_backward and xiMin_forward. (usually 0.1 < epsilon <= 1)
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
xiMin_f                  = 0;
xiMin_b                  = 0;
gamma_Np1_1              = 0;
gamma_N_2                = 0;
gamma_N_3                = 0;
w_f                      = zeros(nCoefficients, 1);
w_b                      = zeros(nCoefficients, 1);
coefficientVector        = zeros(nCoefficients, nIterations+1); 
error_f                  = 0;
error_f_line             = 0;
error_b_line_1           = 0;
error_b_line_2           = 0;
error_b_line_3_Vector    = zeros(3,1);
error_b_3_Vector         = zeros(2,1);
posterioriErrorVector    = zeros(nIterations,1);
prioriErrorVector        = zeros(nIterations,1);
phiHatN                  = zeros(nCoefficients,1);  
phiHatNp1                = zeros(nCoefficients+1,1);
regressor                = zeros(nCoefficients+1,1);


%################################################
% Computations
%################################################

%%Initialize Parameters
w_f                      = zeros(nCoefficients, 1);
w_b                      = zeros(nCoefficients, 1);
coefficientVector(:,1)   = zeros(nCoefficients, 1);
%
phiHatN                  = zeros(nCoefficients,1);  
gamma_N_3                = 1;
%
xiMin_f                  = S.epsilon;
xiMin_b                  = S.epsilon;
%
kappa1                  = 1.5;
kappa2                  = 2.5;
kappa3                  = 1;

prefixedInput            = [ zeros(1,nCoefficients)  input];

%-------------------------------------------------------------
%||||||||||||||||||||||| - Main Loop - |||||||||||||||||||||||
%-------------------------------------------------------------

for it = 1:nIterations
    
    regressor    = prefixedInput(it+nCoefficients:-1:it).';
    
    error_f_line = regressor.'*[ 1 ;  -w_f];
    error_f      = error_f_line*gamma_N_3;
    phiHatNp1    = [0 ; phiHatN ]  + ...
                   1/(S.lambda*xiMin_f)*[1 ; -w_f]*error_f_line;
    
    % Forward Info 
    gamma_Np1_1  = 1/(1/gamma_N_3  +phiHatNp1(1)*error_f_line );
    xiMin_f      = 1/(1/(xiMin_f*S.lambda) -gamma_Np1_1*(phiHatNp1(1))^2);
    w_f          = w_f  +  phiHatN*error_f;
    
    % Backward Errors
    error_b_line_1           = S.lambda*xiMin_b*phiHatNp1(end);
    error_b_line_2           = [-w_b.' 1]*regressor;
    error_b_line_3_Vector(1) = error_b_line_2*kappa1 + ...
                               error_b_line_1*(1 -kappa1);
    error_b_line_3_Vector(2) = error_b_line_2*kappa2 + ...
                               error_b_line_1*(1 -kappa2);
    error_b_line_3_Vector(3) = error_b_line_2*kappa3 + ...
                               error_b_line_1*(1 -kappa3);
    
    % Backward Coefficients
    gamma_N_2 = 1/(1/gamma_Np1_1 - ...
                   phiHatNp1(end)*error_b_line_3_Vector(3));                         
    error_b_3_Vector(1)      =  error_b_line_3_Vector(1)*gamma_N_2;         
    error_b_3_Vector(2)      =  error_b_line_3_Vector(2)*gamma_N_2;
    xiMin_b   =  S.lambda*xiMin_b + ...
                 error_b_3_Vector(2)*error_b_line_3_Vector(2);
    phiHatN   =  phiHatNp1(1:end-1) +phiHatNp1(end)*w_b;
    w_b       =  w_b   +  phiHatN*error_b_3_Vector(1);        
    gamma_N_3 =  1/(1  + phiHatN.'*regressor(1:end-1));
            
    % Joint Process Estimation
    prioriErrorVector(it)     = desired(it) -...
                                coefficientVector(:,it).'*regressor(1:end-1);     
    posterioriErrorVector(it) = prioriErrorVector(it)*gamma_N_3;                             
    coefficientVector(:,it+1) = coefficientVector(:,it) + ...
                                phiHatN*posterioriErrorVector(it);
    

                                
end % END OF LOOP



%EOF

