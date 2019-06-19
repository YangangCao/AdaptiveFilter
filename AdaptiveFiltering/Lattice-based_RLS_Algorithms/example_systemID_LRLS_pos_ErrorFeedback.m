%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        Example: System Identification                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                               %
%  In this example we have a typical system identification scenario. We want    %
% to estimate the filter coefficients of an unknown system given by Wo. In      %
% order to accomplish this task we use an adaptive filter with the same         %
% number of coefficients, N, as the unkown system. The procedure is:            %
% 1)  Excitate both filters (the unknown and the adaptive) with the signal      %
%   x. In this case, x is chosen according to the BPSK constellation.           %
%   The variance of x is 1.                                                     %
% 2)  Generate the desired signal, d = Wo' x + n, which is the output of the    %
%   unknown system considering some disturbance (noise) in the model. The       %
%   noise power is given by sigma_n2.                                           %
% 3)  Choose an adaptive filtering algorithm to govern the rules of coefficient %
%   updating.                                                                   %
%                                                                               %
%     Adaptive Algorithm used here: LRLS_pos_ErrorFeedback                      %
%                                                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%   Definitions:
ensemble    = 100;                          % number of realizations within the ensemble
K           = 500;                          % number of iterations
H           = [0.32+0.21*j,-0.3+0.7*j,0.5-0.8*j,0.2+0.5*j].';
Wo          = real(H);                      % unknown system
sigma_n2    = 0.04;                         % noise power
lambda      = 0.97;                         % forgetting factor
N           = 4;                            % number of coefficients of the adaptive filter (= nSectionsLattice + 1)
epsilon     = 1e-2;                         % small positive constant


%   Initializing & Allocating memory:
ladder  = zeros(N  ,K+1,ensemble);    % ladder coefficients of the algorithm
kappa   = zeros(N  ,K,ensemble);    % reflection coefficients of the lattice algorithm
e       = zeros(N+1,K,ensemble);    % error matrix   
MSE     = zeros(K,ensemble);        % MSE for each realization
MSEmin  = zeros(K,ensemble);        % MSE_min for each realization


%   Computing:
for l=1:ensemble,

    X       = zeros(N,1);               % input at a certain iteration (tapped delay line)
    d       = zeros(1,K);               % desired signal

    x        = sign(randn(K,1));                 % Creating the input signal - BPSK
    sigma_x2 = var(x);                           % signal power = 1
    n        = sqrt(sigma_n2)*(randn(K,1));      % real noise

    for k=1:K,

        X       =   [x(k,1)
                     X(1:(N-1),1)];              % input signal (tapped delay line)

        d(k)    =   (Wo'*X(:,1))+n(k);           % desired signal

    end

    S   =   struct('lambda',lambda,'nSectionsLattice',(N-1),'epsilon',epsilon);
    [ladder(:,:,l),kappa(:,:,l),e(:,:,l)]  =   LRLS_pos_ErrorFeedback(d,transpose(x),S);                   %%% Algorithm 7.3

    MSE(:,l)    =   MSE(:,l)+( (abs(e(N+1,:,l))).^2 ).';
    MSEmin(:,l) =   MSEmin(:,l)+(abs(n(:))).^2;

end


%   Averaging:
ladder_av = sum(ladder,3)/ensemble;
kappa_av  = sum(kappa,3)/ensemble;
MSE_av    = sum(MSE,2)/ensemble;
MSEmin_av = sum(MSEmin,2)/ensemble;


%   Plotting:
figure,
plot(1:K,10*log10(MSE_av),'-k');
title('Learning Curve for MSE (a posteriori error)');
xlabel('Number of iterations, k'); ylabel('MSE [dB]');

figure,
plot(1:K,10*log10(MSEmin_av),'-k');
title('Learning Curve for MSEmin');
xlabel('Number of iterations, k'); ylabel('MSEmin [dB]');

num = zeros(N+1,K);                % estimate of Wo
for k=1:K
  num(:,k) = latc2tf(-kappa_av(:,k),ladder_av(:,k));
end
figure,
plot(real(num(1,:))),...
title('Evolution of the 1st coefficient');
xlabel('Number of iterations, k'); ylabel('Coefficient');

