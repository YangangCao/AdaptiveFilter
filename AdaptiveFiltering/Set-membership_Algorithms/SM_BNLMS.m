function    [outputVector,...
             errorVector,...
             coefficientVector, nUpdates] =   SM_BNLMS(desired,input,S)

%   SM_BNLMS.m
%       Implements the Set-membership Binormalized LMS algorithm for REAL valued data.
%       (Algorithm 6.5 - book: Adaptive Filtering: Algorithms and Practical
%                                                       Implementation, Diniz)
%
%   Syntax:
%       [outputVector,errorVector,coefficientVector,nUpdates] = SM_BNLMS(desired,input,S)
%
%   Input Arguments:
%       . desired   : Desired signal.                               (ROW vector)
%       . input     : Signal fed into the adaptive filter.          (ROW vector)
%       . S         : Structure with the following fields
%           - gamma_bar             : Upper bound for the error modulus.
%           - filterOrderNo         : Order of the FIR filter.
%           - initialCoefficients   : Initial filter coefficients.  (COLUMN vector)
%           - gamma                 : Regularization factor.
%                                     (small positive constant to avoid singularity)
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

%   Pre Allocations
errorVector             =   zeros(nIterations   ,1);
outputVector            =   zeros(nIterations   ,1);
coefficientVector       =   zeros(nCoefficients ,(nIterations+1));
regressor               =   zeros(nCoefficients ,1);

%   Initial State Weight Vector
coefficientVector(:,1)  =   S.initialCoefficients;

%   Improve source code regularity
prefixedInput           =   [zeros(nCoefficients-1,1)
                             transpose(input)];

%   Body
for it = 1:nIterations,

    regressor_prev              =   regressor;

    regressor                   =   prefixedInput(it+(nCoefficients-1):-1:it,1);

    outputVector(it,1)          =   (coefficientVector(:,it).')*regressor;

    errorVector(it,1)           =   desired(it)-outputVector(it,1);

    if abs(errorVector(it,1)) > S.gamma_bar

       mu = 1-(S.gamma_bar/abs(errorVector(it,1)));
       nUpdates = nUpdates+1;

    else

       mu = 0;

    end

    lambda1 = (mu*errorVector(it,1)*(norm(regressor_prev,2))^2)/...
              (S.gamma+((norm(regressor,2))^2)*((norm(regressor_prev,2))^2)-...
              (regressor_prev.'*regressor)^2);

    lambda2 =-(mu*errorVector(it,1)*(regressor_prev.'*regressor))/...
              (S.gamma+((norm(regressor,2))^2)*((norm(regressor_prev,2))^2)-...
              (regressor_prev.'*regressor)^2);


    coefficientVector(:,it+1)   =   coefficientVector(:,it)+...
                                    lambda1*regressor+...
                                    lambda2*regressor_prev;

end

%   EOF
