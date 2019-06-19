function [outputVector, ...
          errorVector,...
          thetaVector,...
          errorVector_e] = ErrorEquation(desired,input,S)

%   ErrorEquation.m
%       Implements the Error Equation RLS algorithm for REAL valued data.
%       (Algorithm 10.3 - book: Adaptive Filtering: Algorithms and Practical
%                                                        Implementation, 3rd Ed., Diniz)
% 
%   Syntax:
%       [outputVector,errorVector,thetaVector] = ErrorEquation(desired,input,S)
% 
%   Input Arguments: 
%       . desired           : Desired Signal.                          (ROW Vector)
%       . input             : Signal fed into the adaptive filter.     (ROW Vector)
%       . S                 : Structure with the following fields
%           - lambda             : Forgetting factor.                  (0 << lambda < 1) 
%           - M                  : Adaptive filter numerator order, refered as M in the textbook.
%           - N                  : Adaptive filter denominator order, refered as N in the textbook.
%           - delta              : Regularization factor. 
% 
%   Output Arguments:
%       . outputVector      : Store the estimated output for each iteration.           (COLUMN vector)
%       . errorVector       : Store the error for each iteration.                      (COLUMN vector)
%       . thetaVector       : Store the estimated coefficients of the IIR filter for each iteration.
%                             (Coefficients at one iteration are COLUMN vector)
%       . errorVector_e     : Store the auxiliary error used for updating thetaVector  (COLUMN vector)
%
%   Authors:
%       . Guilherme de Oliveira Pinto   - guilhermepinto7@gmail.com  &  guilherme@lps.ufrj.br
%       . Markus Vinícius Santos Lima   - mvsl20@gmailcom            &  markus@lps.ufrj.br
%       . Wallace Alves Martins         - wallace.wam@gmail.com      &  wallace@lps.ufrj.br
%       . Luiz Wagner Pereira Biscainho - cpneqs@gmail.com           &  wagner@lps.ufrj.br
%       . Paulo Sergio Ramirez Diniz    -                               diniz@lps.ufrj.br
%

        

% Initialization Procedure
nCoefficients           = S.M + 1 + S.N;  
nIterations             = length(desired);

% Pre Allocations
errorVector             = zeros(nIterations  ,1);
errorVector_e           = zeros(nIterations  ,1);
outputVector            = zeros(nIterations  ,1);
outputVector_e          = zeros(nIterations  ,1);
thetaVector             = zeros(nCoefficients,nIterations+1);
regressor               = zeros(nCoefficients,1);
regressor_e             = zeros(nCoefficients,1);
S_d                     = zeros(nCoefficients,nCoefficients);


% Initial State Weight Vector
%  thetaVector(:,1)         = S.initialCoefficient;
S_d                     = inv(S.delta)*eye(nCoefficients);


% Improve source code regularity (initial state = 0)
prefixedInput           = [zeros(S.M,1)
                           transpose(input)    ]; 
prefixedOutput          = [zeros(S.N,1)
                           outputVector        ]; 
prefixedDesired         = [zeros(S.N,1)
                           transpose(desired)  ]; 


for it = 1:nIterations
    
   regressor                = [ prefixedOutput(it+(S.N-1):-1:it)
                                prefixedInput(it+S.M:-1:it)];
   regressor_e              = [ prefixedDesired(it+(S.N-1):-1:it)
                                prefixedInput(it+S.M:-1:it)];
                              
   % Compute Output                         
   prefixedOutput(it+S.N)   = thetaVector(:,it).'*regressor;
   prefixedOutput_e(it+S.N) = thetaVector(:,it).'*regressor_e;

   % Error
   errorVector(it)          = desired(it) - prefixedOutput(it+S.N);
   errorVector_e(it)        = desired(it) - prefixedOutput_e(it+S.N);        
   
   
   % Sd
   S_d       = inv(S.lambda)*...
               (S_d -(S_d*regressor_e*regressor_e.'*S_d)/...
               (S.lambda + regressor_e.'*S_d*regressor_e));
            
   
    % Update Coefficients
    thetaVector(:,it+1) = thetaVector(:,it) + S_d*regressor_e*errorVector_e(it);                      
    
                                
    %Stability Procedure
    thetaVector(1:S.N, it+1) = stabilityProcedure(thetaVector(1:S.N,it+1));
  
                  
end



outputVector = prefixedOutput(S.N+1:end);



%
% Sub function
%
%


function [coeffOut] = stabilityProcedure(coeffIn)
   

poles                = roots([1 -coeffIn.']);

indexVector          = find(abs(poles) > 1);
%poles(indexVector)   = poles(indexVector)./(abs(poles(indexVector)).^2);
poles(indexVector)   = 1./poles(indexVector);
coeffOut             = poly(poles);

coeffOut              = -real(coeffOut(2:end)).';


%EOF


