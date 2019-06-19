function    [outputVector,...
             errorVector,...
             coefficientVector,...
             outputVectorPost,...
             errorVectorPost] =   RLS(desired,input,S)

%   RLS.m
%       Implements the RLS algorithm for COMPLEX valued data.
%       (Algorithm 5.3 - book: Adaptive Filtering: Algorithms and Practical
%                                                       Implementation, Diniz)
%
%   Syntax:
%       [outputVector,errorVector,coefficientVector,outputVectorPost,...
%                                             errorVectorPost] = RLS(desired,input,S)
%
%   Input Arguments:
%       . desired   : Desired signal.                               (ROW vector)
%       . input     : Signal fed into the adaptive filter.          (ROW vector)
%       . S         : Structure with the following fields
%           - filterOrderNo         : Order of the FIR filter.
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
%   Comments:
%       The authors suggest the use of RLS_Alt function instead of this one. The reason is that
%       the former performs less computations than the latter. Moreover RLS_Alt allows the user
%       to give its own initial coefficient vector.
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
%       . p_d                   :   Cross correlation vector.
%                                   (COLUMN vector)


%   Initialization Procedure
nCoefficients           =   S.filterOrderNo+1;
nIterations             =   length(desired);

%   Pre Allocations
errorVector             =   zeros(nIterations   ,1);
outputVector            =   zeros(nIterations   ,1);
coefficientVector       =   zeros(nCoefficients ,(nIterations+1));

%   Initial State
S_d                     =   S.delta*eye(nCoefficients);
p_d                     =   zeros(nCoefficients,1);

%   Improve source code regularity
prefixedInput           =   [zeros(nCoefficients-1,1)
                             transpose(input)];

%   Body
coefficientVector(:,1)    =   S_d*p_d;

for it = 1:nIterations,

    regressor    =  prefixedInput(it+(nCoefficients-1):-1:it);

    %   a priori estimated output
    outputVector(it,1)  =   coefficientVector(:,it)'*regressor;

    %   a priori error
    errorVector(it,1)   =   desired(it)-outputVector(it,1);

    S_d          =  inv(S.lambda)*(S_d-(S_d*regressor*regressor'*S_d)/...
                    (S.lambda+regressor'*S_d*regressor));

    p_d          =  S.lambda*p_d+conj(desired(it))*regressor;

    coefficientVector(:,it+1)   =   S_d*p_d;

    %    A posteriori estimated output
    outputVectorPost(it,1)      =   coefficientVector(:,it+1)'*regressor;

    %    A posteriori error
    errorVectorPost(it,1)       =   desired(it)-outputVectorPost(it,1);

end

%   EOF
