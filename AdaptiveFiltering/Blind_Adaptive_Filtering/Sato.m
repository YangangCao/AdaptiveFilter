function    [outputVector,...
             errorVector,...
             coefficientVector] =   Sato(input,S)

%   Sato.m
%       Implements the Sato algorithm for COMPLEX valued data.
%       (Algorithm 13.3 - book: Adaptive Filtering: Algorithms and Practical
%                                                        Implementation, Diniz)
%
%   Syntax:
%       [outputVector,errorVector,coefficientVector] = Sato(input,S)
%
%   Input Arguments:
%       . input     : Signal fed into the adaptive filter.          (ROW vector)
%       . S         : Structure with the following fields
%           - step                  : Convergence (relaxation) factor.
%           - filterOrderNo         : Order of the FIR filter.
%           - initialCoefficients   : Initial filter coefficients.  (COLUMN vector)
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
%       . prefixedInput         :   Input is prefixed by nCoefficients -1 random values.
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
%       . desiredLevel          :   Defines the level which abs(outputVector(it,1))^2 
%                                   should approach.


%   Initialization Procedure
nCoefficients       =   S.filterOrderNo+1;
nIterations         =   length(input);
desiredLevel        =   mean(abs(input).^2)/mean(abs(input));

%   Pre Allocations
errorVector             =   zeros(nIterations   ,1);
outputVector            =   zeros(nIterations   ,1);
coefficientVector       =   zeros(nCoefficients ,(nIterations+1));

%   Initial State Weight Vector
coefficientVector(:,1)  =   S.initialCoefficients;

%   Improve source code regularity
prefixedInput           =   [randn(nCoefficients-1,1)
                             transpose(input)];

%   Body
for it = 1:nIterations,

    regressor                   =   prefixedInput(it+(nCoefficients-1):-1:it,1);

    outputVector(it,1)          =   (coefficientVector(:,it)')*regressor;

    errorVector(it,1)           =   outputVector(it,1) - sign(outputVector(it,1))*desiredLevel;

    coefficientVector(:,it+1)   =   coefficientVector(:,it)-...
                                    (S.step*conj(errorVector(it,1))*...
				    regressor);

end

%   EOF
