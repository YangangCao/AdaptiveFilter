function    [outputVector,...
             errorVector,...
             coefficientVector, nUpdates] =   Simp_SM_PUAP(desired,input,S)

%   Simp_SM_PUAP.m
%       Implements the Simplified Set-membership Partial-Update Affine-Projection
%                                               algorithm for COMPLEX valued data.
%       (Algorithm 6.6 - book: Adaptive Filtering: Algorithms and Practical
%                                                       Implementation, Diniz)
%
%   Syntax:
%       [outputVector,errorVector,coefficientVector,nUpdates] = Simp_SM_PUAP(desired,input,S)
%
%   Input Arguments:
%       . desired   : Desired signal.                               (ROW vector)
%       . input     : Signal fed into the adaptive filter.          (ROW vector)
%       . S         : Structure with the following fields
%           - gamma_bar             : Upper bound for the error modulus.
%           - filterOrderNo         : Order of the FIR filter.
%           - initialCoefficients   : Initial filter coefficients.              (COLUMN vector)
%           - upSelector            : Indicates which coefficients will be updated.    (vector)
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
%         The Partial-Update means that when the (adaptive) filter coefficients have to be
%       updated, only some coefficients, given by S.upSelector, are actually changed. So,
%       each row of S.upSelector (for a given iteration, i.e., column) must have only 0 or
%       or 1 elements.
%       Ex:
%                S.upSelector(:,10) = [1 1 0 0].'
%
%       means that if it is necessay to update in the 10th iteration, only the first and
%       the second coefficients are going to be changed. Notice that S.upSelector must be
%       a matrix with (S.filterOrderNo+1) rows and (nIterations) columns.
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
size_upSelector     =   size(S.upSelector);
if size_upSelector(1) ~= nCoefficients
   fprintf('S.upSelector must have N (number of filter coefficients) rows. \n \n');
   return
end
if size_upSelector(2) ~= nIterations
   fprintf('S.upSelector must have K (number of iterations) columns. \n \n');
   return
end

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

    C                          =  diag(S.upSelector(:,it));

    coefficientVector(:,it+1)  =  coefficientVector(:,it)+ C*(regressor*...
                                  inv(regressor'*C*regressor+S.gamma*eye(S.memoryLength+1))*...
                                  mu*errorVectorApConj(1,it)*u1);

end

outputVector            =   outputVectorApConj(1,:)';
errorVector             =   errorVectorApConj(1,:)';

%   EOF
