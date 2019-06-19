%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adaptive Filter Theory 5e Solution Manual                       %
%                                                                 %
% Chapter 6                                                       %
% Question 21                                                     %
%                                                                 %
% Program written to run on MATLAB 2010a (R)                      %
%                                                                 %
% By Kelvin Hall                                                  %
% July 2, 2014                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc
numberOfDatapoints=1500; % try 300 to get a better picture of the process
numberOfRuns=200;       % increase this value to get smoothed out image

W=[2.9,3.1,3.5]; %eigen value spread parameter
NoiseVariance=0.001; % the variance of the system noise given in the problem
mu=[0.0075, 0.025, 0.075]; %LMS learning rate parameters

stream = RandStream('mt19937ar','Seed',3);  % seed the random number
RandStream.setDefaultStream(stream);        % generator for reproducable
                                            % results
NoiseStandardDeviation=sqrt(NoiseVariance); 

u=zeros(numberOfDatapoints+14,1);  % Allocate memory for input data stream
r=zeros(numberOfDatapoints+14,1);  % Allocate memory for recieved data
f=zeros(numberOfDatapoints+14,1);  % Allocate memory for error between 
                                % prediction and resutls
g=zeros(numberOfDatapoints+14,1);  % allocate memory of squared averaged error
for ell=1:3
    h=[0.5*(1+cos(2*pi/W(ell)*(-1))),0.5*(1+cos(2*pi/W(ell)*(0))),...
                                        0.5*(1+cos(2*pi/W(ell)*(1)))];
    for m=1:size(mu,2)
    for k=1:numberOfRuns % Loop for performing the appropriate number of 
                     % Monte Carlo simulations
    u=zeros(numberOfDatapoints,1);          % Initialize the weights to zero
    
    for n=3+11:numberOfDatapoints+3+11 % The loop that does the data runs and 
                                %filter updates
 
        u(n)=binornd(1,0.5)*2-1;
        r(n)=u(n)*h(3)+u(n-1)*h(2)+u(n-2)*h(1)+randn(1)*NoiseStandardDeviation;
        f(n)=u(n-5)-Weight'*r(n:-1:n-10); % calculate the error in estimation
        Weight=Weight+mu(m)*r(n:-1:n-10)*f(n); % update the weights in the manner associated
                            % LMS based on the recient error in prediction
                            % of the last input
    end
    g=g+f.^2;               % accumulate squared error of estimation
    end
    g=g/numberOfRuns;           % normalize the accumulated error to reflect
                            % the average error based on the number of
                            % monte Carlo simulations completed
    subplot(3,3,(ell-1)*3+m)
    semilogy([1:numberOfDatapoints+14],g,'k') % plot the MSE vs iterations
    str=sprintf('Graph of squared error vs iteration \n with mu= %.4f and W=%.1f'...
                ,mu(m),W(ell));
    title(str)
    xlabel('Number of iterations') % x-axis label
    xlim([0 1500])
    ylim([1e-3 1])
    ylabel('Squared Error') % y-axis label
    end
end