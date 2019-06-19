function    [outputVector,...
             errorVector,...
             coefficientVector] =   Dual_sign(desired,input,S)

%   Dual_sign.m
%       Implements the Dual-Sign LMS algorithm for REAL valued data.
%       (Modified version of Algorithm 4.1 - book: Adaptive Filtering: Algorithms
%                                               and Practical Implementation, Diniz)
%
%   Syntax:
%       [outputVector,errorVector,coefficientVector] = Dual_sign(desired,input,S)
%
%   Input Arguments:
%       . desired   : Desired signal.                               (ROW vector)
%       . input     : Signal fed into the adaptive filter.          (ROW vector)
%       . S         : Structure with the following fields
%           - step                  : Convergence (relaxation) factor.(Factor 2 included)
%           - filterOrderNo         : Order of the FIR filter.
%           - initialCoefficients   : Initial filter coefficients.  (COLUMN vector)
%           - rho                   : Error modulus threshold.
%           - gamma                 : Gain factor. Gamma is a power of two. (gamma > 1)
%
%   Output Arguments:
%       . outputVector      :   Store the estimated output of each iteration.   (COLUMN vector)
%       . errorVector       :   Store the error for each iteration.             (COLUMN vector)
%       . coefficientVector :   Store the estimated coefficients for each iteration.
%                               (Coefficients at one iteration are COLUMN vector)
%
%   Authors:
%       . Guilherme de Oliveira Pinto   - guilhermepinto7@gmail.com & guilherme@lps.ufrj.br
%       . Markus Vinícius Santos Lima   - mvsl20@gmailcom           & markus@lps.ufrj.br
%       . Wallace Alves Martins         - wallace.wam@gmail.com     & wallace@lps.ufrj.br
%       . Luiz Wagner Pereira Biscainho - cpneqs@gmail.com          & wagner@lps.ufrj.br
%       . Paulo Sergio Ramirez Diniz    -                             diniz@lps.ufrj.br
%


%   Some Variables and Definitions:
%       . prefixedInput         :   Input is prefixed by nCoefficients -1 zeros.
%                                   (The prefix led to a more regular source code)
%
%       . regressor             :   Auxiliar variable. Store the piece of the
%                                   prefixedInput that will be multiplied by the
%                                   current set of coefficients.
%                                   (regressor is a COLUMN vector)
%
%       . nCoefficients         :   FIR filter number of coefficients.
%
%       . nIterations           :   Number of iterations.


%   Initialization Procedure
nCoefficients       =   S.filterOrderNo+1;
nIterations         =   length(desired);

%   Pre Allocations
errorVector             =   zeros(nIterations   ,1);
outputVector            =   zeros(nIterations   ,1);
coefficientVector       =   zeros(nCoefficients ,(nIterations+1));
dualSignError           =   0;

%   Initial State Weight Vector
coefficientVector(:,1)  =   S.initialCoefficients;

%   Improve source code regularity
prefixedInput           =   [zeros(nCoefficients-1,1)
                             transpose(input)];

%   Body
for it=1:nIterations,

    regressor                   =   prefixedInput(it+(nCoefficients-1):-1:it,1);

    outputVector(it,1)          =   (regressor.')*coefficientVector(:,it);

    errorVector(it,1)           =   desired(it)-outputVector(it,1);

    if(abs(errorVector(it,1))>S.rho)
        dualSignError           =   S.gamma*sign(errorVector(it,1));
    else
        dualSignError           =   sign(errorVector(it,1));
    end

    coefficientVector(:,it+1)   =   coefficientVector(:,it)+(...
                                    2*S.step*dualSignError*...
                                    regressor);

end

%   EOF
