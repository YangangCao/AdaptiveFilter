%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        Example: System Identification                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                               %
%  In this example we have a typical system identification scenario. We want    %
% to estimate the filter coefficients of an unknown system given by Wo. In      %
% order to accomplish this task we use an adaptive filter with the same         %
% number of coefficients, N, as the unkown system. The procedure is:            %
% 1)  Excitate both filters (the unknown and the adaptive) with the signal      %
%   x. In this case, x is chosen according to the 4-QAM constellation.          %
%   The variance of x is normalized to 1.                                       %
% 2)  Generate the desired signal, d = Wo' x + n, which is the output of the    %
%   unknown system considering some disturbance (noise) in the model. The       %
%   noise power is given by sigma_n2.                                           %
% 3)  Choose an adaptive filtering algorithm to govern the rules of coefficient %
%   updating.                                                                   %
%                                                                               %
%     Adaptive Algorithm used here: ErrorEquation                               %
%                                                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%   Definitions:
ensemble    = 100;                          % number of realizations within the ensemble
K           = 500;                          % number of iterations
H           = [0.32+0.21*j,-0.3+0.7*j,0.5-0.8*j,0.2+0.5*j].';
Wo          = real(H);                      % unknown system
sigma_n2    = 0.001;                         % noise power
M           = 3;                            % number of coefficients of the adaptive filter
N           = 2;
lambda      = 0.97;                         % forgetting factor
delta       = 1e-5;


%   Initializing & Allocating memory:
theta   = zeros(M+1+N,K+1,ensemble);  % coefficient vector for each iteration and realization
MSE     = zeros(K,ensemble);          % MSE for each realization
MSEmin  = zeros(K,ensemble);          % MSE_min for each realization


%   Computing:
for l=1:ensemble,

    X       = zeros(M+1,1);               % input at a certain iteration (tapped delay line)
    d       = zeros(1,K);                 % desired signal

    x        = sign(randn(K,1));          % Creating the input signal (normalized)
    sigma_x2 = var(x);                    % signal power = 1

    n        = sqrt(sigma_n2)*randn(K,1);       % complex noise

    for k=1:K,

        X       =   [x(k,1)
                     X(1:(M),1)];              % input signal (tapped delay line)

        d(k)    =   (Wo'*X(:,1))+n(k);           % desired signal

    end

    S   =   struct('lambda',lambda,'M',M,'N',N,'delta',delta);
    [y,e,theta,ee]  =   ErrorEquation(d,transpose(x),S);


    MSE(:,l)    =   MSE(:,l)+(abs(e(:,1))).^2;
    MSEE(:,l)   =   MSE(:,l)+(abs(ee(:,1))).^2;
    MSEmin(:,l) =   MSEmin(:,l)+(abs(n(:))).^2;

end


%   Averaging:
theta_av    = sum(theta,3)/ensemble;
MSE_av      = sum(MSE,2)/ensemble;
MSEE_av     = sum(MSEE,2)/ensemble;
MSEmin_av   = sum(MSEmin,2)/ensemble;


%   Plotting:
figure,
plot(1:K,10*log10(MSE_av),'-k');
title('Learning Curve for MSE');
xlabel('Number of iterations, k'); ylabel('MSE [dB]');

figure,
plot(1:K,10*log10(MSEmin_av),'-k');
title('Learning Curve for MSEmin');
xlabel('Number of iterations, k'); ylabel('MSEmin [dB]');


fprintf('Adaptive Filter Coefficients (last iteration computed over %d runs): \n',ensemble);
fprintf('Numerator coefficients (direct part): \n ');
theta_av(N+1:end,end).'
fprintf('Denominator coefficients (recursive part): \n ');
theta_av(1:N,end).'
fprintf('Unkown system (FIR): \n ');
Wo.'

