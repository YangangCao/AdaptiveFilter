function [outputVector, ...
          errorVector,...
          thetaVector] = GaussNewton_GradientBased(desired,input,S)

%   GaussNewton_GradientBased.m
%       Implements the Gradient-Based Gauss-Newton algorithm for REAL valued data.
%       (Modified version of Algorithm 10.1 - book: Adaptive Filtering: Algorithms and 
%                                                          Practical Implementation, 3rd Ed., Diniz)
% 
%   Syntax:
%       [outputVector,errorVector,thetaVector] = GaussNewton_GradientBased(desired,input,S)
% 
%   Input Arguments: 
%       . desired           : Desired Signal.                          (ROW Vector)
%       . input             : Signal fed into the adaptive filter.     (ROW Vector)
%       . S                 : Structure with the following fields
%           - alpha              : scalar factor similar to (1 - Forgetting factor).              
%           - step               : Step-size.
%           - M                  : Adaptive filter numerator order, refered as M in the textbook.
%           - N                  : Adaptive filter denominator order, refered as N in the textbook.
% 
%   Output Arguments:
%       . outputVector      : Store the estimated output for each iteration.  (COLUMN vector)
%       . errorVector       : Store the error for each iteration.             (COLUMN vector)
%       . thetaVector       : Store the estimated coefficients of the IIR filter for each iteration.
%                             (Coefficients at one iteration are COLUMN vector)
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
outputVector            = zeros(nIterations  ,1);
thetaVector             = zeros(nCoefficients,nIterations+1);
regressor               = zeros(nCoefficients,1);


% Initial State Weight Vector
%  thetaVector(:,1)         = S.initialCoefficient;
xLineVector             = zeros(max(S.M+1, S.N),   1);
yLineVector             = zeros(max(S.M+1, S.N),   1);
PhiVector               = zeros(nCoefficients,1);


% Improve source code regularity (initial state = 0)
prefixedInput           = [zeros(S.M,1)
                           transpose(input)    ]; 
prefixedOutput          = [zeros(S.N,1)
                           outputVector        ]; 


for it = 1:nIterations
    
   regressor                = [ prefixedOutput(it+(S.N-1):-1:it)
                                prefixedInput(it+S.M:-1:it)];
                              
   % Compute Output                         
   prefixedOutput(it+S.N)   = thetaVector(:,it).'*regressor;

   % Error
   errorVector(it)          = desired(it) - prefixedOutput(it+S.N);    
   
   %Xline,Yline
   xLineK      =   prefixedInput(it+S.M)  + thetaVector(1:S.N,it).'*xLineVector(1:S.N); 
   yLineK      =  -prefixedOutput(it+S.N-1) + thetaVector(1:S.N,it).'*yLineVector(1:S.N); 
   xLineVector = [ xLineK ; xLineVector(1:end-1) ];
   yLineVector = [ yLineK ; yLineVector(1:end-1) ];
   PhiVector   = [ +yLineVector(1:S.N) ; -xLineVector(1:S.M+1) ];
         
      
    % Update Coefficients
    thetaVector(:,it+1) = thetaVector(:,it) -S.step*PhiVector*errorVector(it);                      
    
                                
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


