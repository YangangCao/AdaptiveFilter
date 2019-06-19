function    [outputVector,...
             errorVector,...
             coefficientVector, nUpdates] =   Simp_SM_AP(desired,input,S)

%   Simp_SM_AP.m
%       Implements the Simplified Set-membership Affine-Projection algorithm for
%                                                            COMPLEX valued data.
%       (Algorithm 6.7 - book: Adaptive Filtering: Algorithms and Practical
%                                                       Implementation, Diniz)
%
%   Syntax:
%       [outputVector,errorVector,coefficientVector,nUpdates] = Simp_SM_AP(desired,input,S)
%
%   Input Arguments:
%       . desired   : Desired signal.                               (ROW vector)
%       . input     : Signal fed into the adaptive filter.          (ROW vector)
%       . S         : Structure with the following fields
%           - gamma_bar             : Upper bound for the error modulus.
%           - filterOrderNo         : Order of the FIR filter.
%           - initialCoefficients   : Initial filter coefficients.              (COLUMN vector)
%           - gamma                 : Regularization factor.
%                                     (small positive constant to avoid singularity)
%           - memoryLength          : Reuse data factor.
%                                     (referred as L in the textbook)
%
%   Output Arguments:
%       . outputVector      :   Store the estimated output of each iteration.   (COLUMN vector)
%       . errorVector       :   Store the error for each iteration.             (COLUMN vector)
%       . coefficientVector :   Store the estimated coefficients for each iteration.
%                               (Coefficients at one iteration are COLUMN vector)
%       . nUpdates          :   Number of filter coefficient updates.
%
%   Comments:
%         Set-membership filtering implies that the (adaptive) filter coefficients are only
%       updated if the magnitude of the error is greater than S.gamma_bar. In practice, we
%       choose S.gamma_bar as a function of the noise variance (sigma_n2). A commom choice
%       is S.gamma_bar = sqrt(5 * sigma_n2).
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
nUpdates            =   0;
u1                  =   [1 zeros(1,S.memoryLength)].';

%   Pre Allocations
errorVectorApConj       =   zeros((S.memoryLength+1)    ,nIterations);
outputVectorApConj      =   zeros((S.memoryLength+1)    ,nIterations);
coefficientVector       =   zeros(nCoefficients         ,(nIterations+1));
regressor               =   zeros(nCoefficients         ,(S.memoryLength+1));

%   Initial State Weight Vector
coefficientVector(:,1)  =   S.initialCoefficients;

%   Improve source code regularity
prefixedInput           =   [zeros((nCoefficients-1),1)
                             transpose(input)];
prefixedDesired         =   [zeros(S.memoryLength,1)
                             transpose(desired)];

%   Body
for it = 1:nIterations,

    regressor(:,2:S.memoryLength+1) =   regressor(:,1:S.memoryLength);

    regressor(:,1)                  =   prefixedInput(it+(nCoefficients-1):-1:it);

    outputVectorApConj(:,it)        =   (regressor')*coefficientVector(:,it);

    errorVectorApConj(:,it)         =   conj(prefixedDesired(it+(S.memoryLength):-1:it))...
                                        -outputVectorApConj(:,it);

    if abs(errorVectorApConj(1,it)) > S.gamma_bar

       nUpdates = nUpdates+1;
       mu = 1-(S.gamma_bar/abs(errorVectorApConj(1,it)));

    else

       mu = 0;

    end

    coefficientVector(:,it+1)  =  coefficientVector(:,it)+(regressor*...
                                  inv(regressor'*regressor+S.gamma*eye(S.memoryLength+1))*...
                                  mu*errorVectorApConj(1,it)*u1);

end

outputVector            =   outputVectorApConj(1,:)';
errorVector             =   errorVectorApConj(1,:)';

%   EOF
