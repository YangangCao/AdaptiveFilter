%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adaptive Filter Theory 5e Solution Manual                           %
%                                                                     %
% Chapter 10                                                          %
% Question 10                                                         %
%                                                                     %
%                                                                     %
% Program written to run on MATLAB 2010a (R)                          %
%                                                                     %
% By Kelvin Hall                                                      %
% June 30, 2014                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
lambda=0.9;   % forgetting factor
delta=0.005;  % regularization paramter
H1=[0.25,1,0.25];   %channel response one
H2=[0.25,1,-0.25];  %channel response two
H3=[-0.25,1,0.25];  %channel response three
standardDeviationNoise=sqrt(0.01);  %addivite white noise standard deviation
delay=11;       % optimal delay found in earlier question (ch6 question 18)

numberOfDatapoints=600;     % number of data points per run
NumberOfRuns=200;       % number of monte carlo runs that the results 
                        %are averaged over
u=zeros(numberOfDatapoints+delay+1+21,1);  % allocate memory for filter input
r=zeros(1,numberOfDatapoints+delay+1+21);  % allocate memory for filter output
error=zeros(numberOfDatapoints+delay+1+21,1);  %allocate memory for errors
E=zeros(numberOfDatapoints+delay+1+21,1); % allocate memory for MSE

for k=1:NumberOfRuns
    P=eye(21)*delta; %reset filter variable between monte carlo runs
    W=zeros(21,1);
    for n=delay+1+21:numberOfDatapoints+delay+1+21
        u(n)=rand(1); % create new input data
        if u(n)>0.5
            u(n)=1;
        else
            u(n)=-1;
        end
        d=u(n-delay);  
        r(n)=H1*[u(n);u(n-1);u(n-2)]+randn*standardDeviationNoise;
        
        U=r(n:-1:n-20)'; % create a vector of the 21 most recent filter outputs
        
        error(n)=u(n-delay)-U'*W; 
        
        kappa=lambda^(-1)*P*U/(1+lambda^(-1)*U'*P*U); % update kappa as perRLS
        
        W=W+kappa*error(n); % update weights
        
        P=lambda^(-1)*P-lambda^(-1)*kappa*U'*P; %update as per RLS
    end
    E=E+error.^2; % accumulate squared error in prediction
end
E=E/NumberOfRuns; % divide the squared error by the number of 
                 % monte carlo runs to get the Mean Squared Error

semilogy(E)
ylim([0.01 1])
xlim([45 625])
title('Learning Curve')
xlabel('Iterations')
ylabel('Mean Squared Error')

    