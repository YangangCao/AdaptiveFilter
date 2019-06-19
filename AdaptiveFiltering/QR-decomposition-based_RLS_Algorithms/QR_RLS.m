function [errorVector,...
          coefficientVector] = QR_RLS(desired,input,S)

%   QR_RLS.m
%       Implements the QR-RLS algorithm for REAL valued data.
%       (Algorithm 9.1 - book: Adaptive Filtering: Algorithms and Practical
%                                                       Implementation, 3rd Ed., Diniz)
% 
%   Syntax:
%       [errorVector,coefficientVector] = QR_RLS(desired,input,S)
% 
%   Input Arguments: 
%       . desired           : Desired Signal.                          (ROW Vector)
%       . input             : Signal fed into the adaptive filter.     (ROW Vector)
%       . S                 : Structure with the following fields
%           - filterOrderNo    : Order of the FIR filter.
%           - lambda           : forgetting factor.                    (0 << lambda < 1)
% 
%   Output Arguments:
%       . errorVector       : Store the A POSTERIORI error for each iteration.    (COLUMN Vector)
%       . coefficientVector : Store the estimated coefficients for each iteration.
%                             (Coefficients at one iteration are COLUNN vector)
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
nCoefficients           = S.filterOrderNo +1;  
nIterations             = length(desired); 

%% Pre Allocations
coefficientVector       = zeros(nCoefficients,nIterations);
errorVector             = zeros(nIterations,1);
ULineMatrix             = zeros(nCoefficients,nCoefficients);
dLine_q2                = zeros(nCoefficients,1);
QLineTeta               = zeros(nCoefficients,nCoefficients);
auxUline                = zeros(nCoefficients,nCoefficients); 
auxDline                = zeros(nCoefficients,nCoefficients); 
regressorLine           = zeros(1,nCoefficients);
dBar                    = zeros(nCoefficients,1);
% Scalar Values
gamma                   = 0;
cosTetaI                = 0;
sinTetaI                = 0;
cI                      = 0;
dLine                   = 0;



%################################################
% Computations
%################################################

%%%%% Backsubstituition Procedure
coefficientVector(1,1) = desired(1)/input(1);
for kt = 2:nCoefficients
    
    coefficientVector(1,kt) = desired(1)/input(1);
    for ct = 2:kt
        
        coefficientVector(ct,kt) = ...
        ( -input(2:ct)*coefficientVector(ct-1:-1:1,kt) + desired(ct))/input(1);
    end
end


%%%%% Build Initial Matrices
ULineMatrix           = zeros(nCoefficients,nCoefficients);
for it = 0:S.filterOrderNo
    %Uline
    ULineMatrix(it+1,1:end-it) = ....
    (S.lambda^((it+1)/2))*input(nCoefficients -it:-1:1);
    
    %dLine_q2
    dLine_q2(it+1) = S.lambda^((it+1)/2)*desired(nCoefficients-it);
                        
end

%-------------------------------------------------------------
%||||||||||||||||||||||| - Main Loop - |||||||||||||||||||||||
%-------------------------------------------------------------

for kt = nCoefficients+1:nIterations

    gamma                   = 1;
    dLine                   = desired(kt);
    regressorLine           = input(kt:-1:kt -S.filterOrderNo);
    
    % Givens Rotations 
    for rt = 0:S.filterOrderNo
    
        %Sine and cosine parameters
        cI        = sqrt(ULineMatrix(rt+1,nCoefficients-rt)^2 + ...
                         regressorLine(end-rt)^2 );
        cosTetaI  = ULineMatrix(rt+1,nCoefficients-rt)/cI;
        sinTetaI  = regressorLine(end-rt)/cI;
        
        %Rotation Matrix
        auxInd    = S.filterOrderNo -rt;
 QLineTeta = [cosTetaI        zeros(1,rt)     -sinTetaI        zeros(1,auxInd)
              zeros(rt,1)     eye(rt)          zeros(rt,1)     zeros(rt,auxInd)  
              sinTetaI        zeros(1,rt)      cosTetaI        zeros(1,auxInd)
              zeros(auxInd,1) zeros(auxInd,rt) zeros(auxInd,1) eye(auxInd) ];
                 
        %% Givens Rotations
        %U
        auxUline      = QLineTeta*[regressorLine ; ULineMatrix];            
        regressorLine = auxUline(1,:);
        ULineMatrix   = auxUline(2:end,:);
        gamma = gamma*cosTetaI;
        %d
        auxDline = QLineTeta*[dLine  ; dLine_q2];
        dLine = auxDline(1,:);
        dLine_q2 = auxDline(2:end,:);
        
        
        
    end
   
    %Compute Coeffcient Vector
    dBar = [dLine ; dLine_q2];
    coefficientVector(1,kt) = dBar(nCoefficients+1)/...
                              ULineMatrix(nCoefficients,1);
    for it =1:S.filterOrderNo

        coefficientVector(it+1,kt) = ...
       ( -ULineMatrix(nCoefficients-it,it:-1:1)*coefficientVector(it:-1:1,kt) + dBar(nCoefficients+1-it))/...
       ULineMatrix(nCoefficients-it,1+it);
   
    end

    %Update Parameters 
    dLine_q2        = (S.lambda^0.5)*dLine_q2;
    ULineMatrix     = (S.lambda^0.5)*ULineMatrix;
    errorVector(kt) = dLine*gamma;



end % END OF LOOP

   
%EOF
