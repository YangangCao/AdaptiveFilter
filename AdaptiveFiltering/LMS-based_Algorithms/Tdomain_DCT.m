function    [outputVector,...
             errorVector,...
             coefficientVectorDCT,...
             coefficientVector] =   Tdomain_DCT(desired,input,S)

%   Tdomain_DCT.m
%       Implements the Transform-Domain LMS algorithm, based on the Discrete
%       Cossine Transform (DCT) Matrix, for COMPLEX valued data.
%       (Algorithm 4.4 - book: Adaptive Filtering: Algorithms and Practical
%                                                       Implementation, Diniz)
%
%   Syntax:
%       [outputVector,errorVector,coefficientVectorDCT,coefficientVector] = Tdomain_DCT(desired,input,S)
%
%   Input Arguments:
%       . desired   : Desired signal.                               (ROW vector)
%       . input     : Signal fed into the adaptive filter.          (ROW vector)
%       . S         : Structure with the following fields
%           - step                  : Convergence (relaxation) factor.
%           - filterOrderNo         : Order of the FIR filter.
%           - initialCoefficients   : Initial filter coefficients
%                                     in the ORIGINAL domain.       (COLUMN vector)
%           - gamma                 : Regularization factor.
%                                     (small positive constant to avoid singularity)
%           - alpha                 : Used to estimate eigenvalues of Ru.
%                                     (0 << alpha < 0.1)
%           - initialPower          : Initial power.                (SCALAR)
%
%   Output Arguments:
%       . outputVector          :   Store the estimated output of each iteration.
%                                   (COLUMN vector)
%       . errorVector           :   Store the error for each iteration.
%                                   (COLUMN vector)
%       . coefficientVectorDCT  :   Store the estimated coefficients for each iteration
%                                   in the TRANSFORM domain.
%                                   (Coefficients at one iteration are COLUMN vector)
%       . coefficientVector     :   Store the estimated coefficients for each iteration
%                                   in the ORIGINAL domain.
%                                   (Coefficients at one iteration are COLUMN vector)
%
%   Comments:
%       The adaptive filter is implemented in the Transform-Domain (DCT). Therefore, the first
%       three output variables are calculated in this TRANSFORMED domain. The last output
%       variable, coefficientVector, corresponds to the adaptive filter coefficients in the
%       ORIGINAL domain (the coefficientVector is the Inverse Discrete Cossine Transform
%       aplied to the coefficientVectorDCT) and is only calculated in order to facilitate
%       comparisons, i.e., for implementation purposes just coefficientVectorDCT matters.
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
%       . regressorDCT          :   Auxiliar variable. Store the transformed piece
%                                   of the prefixedInput that will be multiplied
%                                   by the current set of transformed coefficients.
%                                   (regressorDCT is a COLUMN vector)
%
%       . nCoefficients         :   FIR filter number of coefficients.
%
%       . nIterations           :   Number of iterations.
%
%       . coefficientVectorDCT  :   Store the DCT transformed set of weight factors.
%                                   (coefficientVectorDCT is a COLUMN vector)


%   Initialization Procedure
nCoefficients       =   S.filterOrderNo+1;
nIterations         =   length(desired);
T                   =   dctmtx(nCoefficients);

%   Pre Allocations
errorVector             =   zeros(nIterations   ,1);
outputVector            =   zeros(nIterations   ,1);
coefficientVectorDCT    =   zeros(nCoefficients ,(nIterations+1));

%   Initial State
coefficientVectorDCT    =   T*(S.initialCoefficients);
powerVector             =   S.initialPower*ones(nCoefficients,1);

%   Improve source code regularity
prefixedInput           =   [zeros(nCoefficients-1,1)
                             transpose(input)];

%   Body
for it = 1:nIterations,

    regressorDCT        =   T*(prefixedInput(it+(nCoefficients-1):-1:it));

    %   Summing two column vectors
    powerVector         =   S.alpha*(regressorDCT.*conj(regressorDCT))+...
                            (1-S.alpha)*(powerVector);

    outputVector(it,1)  =   (coefficientVectorDCT(:,it)')*regressorDCT;

    errorVector(it,1)   =   desired(it)-outputVector(it,1);

    %   Vectorized
    coefficientVectorDCT(:,it+1)=   coefficientVectorDCT(:,it)+(...
                                    (S.step*conj(errorVector(it,1))*...
                                    regressorDCT)./(S.gamma+powerVector));

end

coefficientVector = T'*(coefficientVectorDCT);

%   EOF
