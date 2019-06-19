function    [outputVector,...
             errorVector,...
             coefficientVector,...
             outputVectorPost,...
             errorVectorPost] =   RLS_Alt(desired,input,S)

%   RLS_Alt.m
%       Implements the Alternative RLS algorithm for COMPLEX valued data. RLS_alt
%       differs from RLS in the number of computations. The RLS_alt function uses
%       an auxiliar variable (psi) in order to reduce the computational burden.
%       (Algorithm 5.4 - book: Adaptive Filtering: Algorithms and Practical
%                                                           Implementation, Diniz)
%
%   Syntax:
%       [outputVector,errorVector,coefficientVector,outputVectorPost,...
%                                             errorVectorPost] = RLS_Alt(desired,input,S)
%
%   Input Arguments:
%       . desired   : Desired signal.                               (ROW vector)
%       . input     : Signal fed into the adaptive filter.          (ROW vector)
%       . S         : Structure with the following fields
%           - filterOrderNo         : Order of the FIR filter.
%           - initialCoefficients   : Initial filter coefficients.  (COLUMN vector)
%           - delta                 : The matrix delta*eye is the initial value of the
%                                     inverse of the deterministic autocorrelation matrix.
%           - lambda                : Forgetting factor.            (0 << lambda < 1)
%
%   Output Arguments:
%       . outputVector          : Store the estimated output of each iteration.
%                                                                   (COLUMN vector)
%       . errorVector           : Store the error for each iteration.
%                                                                   (COLUMN vector)
%       . coefficientVector     : Store the estimated coefficients for each iteration.
%                                 (Coefficients at one iteration are COLUMN vector)
%       . outputVectorPost      : Store the a posteriori estimated output of each iteration.
%                                                                   (COLUMN vector)
%       . errorVectorPost       : Store the a posteriori error for each iteration.
%                                                                   (COLUMN vector)
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
%
%       . S_d                   :   Inverse of the deterministic autocorrelation matrix.
%                                   (It is not an estimate since we are working with
%                                   the least squares approach)
%
%       . psi                   :   Auxiliar variable. Improve computations.
%                                   (COLUMN vector)
%


%   Initialization Procedure
nCoefficients           =   S.filterOrderNo+1;
nIterations             =   length(desired);

%   Pre Allocations
errorVector             =   zeros(nIterations   ,1);
outputVector            =   zeros(nIterations   ,1);
coefficientVector       =   zeros(nCoefficients ,(nIterations+1));

%   Initial State
coefficientVector(:,1)  =   S.initialCoefficients;
S_d                     =   S.delta*eye(nCoefficients);

%   Improve source code regularity
prefixedInput           =   [zeros(nCoefficients-1,1)
                             transpose(input)];

%   Body
for it = 1:nIterations,

    regressor           =   prefixedInput(it+(nCoefficients-1):-1:it);

    %   a priori estimated output
    outputVector(it,1)  =   coefficientVector(:,it)'*regressor;

    %   a priori error
    errorVector(it,1)   =   desired(it)-outputVector(it,1);

    psi                 =   S_d*regressor;

    S_d                 =   inv(S.lambda)*(S_d-(psi*psi')/...
                            (S.lambda + psi'*regressor));

    coefficientVector(:,it+1)   =   coefficientVector(:,it)+...
                                    conj(errorVector(it,1))*S_d*regressor;

    %    A posteriori estimated output
    outputVectorPost(it,1)      =   coefficientVector(:,it+1)'*regressor;

    %    A posteriori error
    errorVectorPost(it,1)       =   desired(it)-outputVectorPost(it,1);

end

%   EOF
