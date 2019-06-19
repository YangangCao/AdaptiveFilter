%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adaptive Filter Theory 5e Solution Manual                            %
%                                                                      %
% Chapter 6                                                            %
% Question 19                                                          %
%                                                                      %
% Program written to run on MATLAB 2010a (R)                           %
%                                                                      %
% By Kelvin Hall                                                       %
% July 2, 2014                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc
numberOfDatapoints=100; % try 300 to get a better picture of the process
numberOfRuns=100;       % increase this value to get smoothed out image

a=-0.99; %AR paramter
NoiseVariance=0.02; % the variance of the system noise given in the problem
mu=0.01; %LMS learning rate parameter it will interest the reader to also
                        %try also 0.5 and 0.75

stream = RandStream('mt19937ar','Seed',6);  % seed the random number
RandStream.setDefaultStream(stream);        % generator for reproducable
                                            % results
NoiseStandardDeviation=sqrt(NoiseVariance); 

u=zeros(numberOfDatapoints,1);  % Allocate memory for input data stream
f=zeros(numberOfDatapoints,1);  % Allocate memory for error between 
                                % prediction and resutls
g=zeros(numberOfDatapoints,1);  % allocate memory of squared averaged error

for k=1:numberOfRuns % Loop for performing the appropriate number of 
                     % Monte Carlo simulations
    u(1)=0;          % Initialize the weights back to zero as well as the
    w=0;             % initial input to zero
    
    for n=2:numberOfDatapoints % The loop that does the data runs and 
                                %filter updates
 
        u(n)=-a*u(n-1)+randn(1)*NoiseStandardDeviation;% update input data 
                            % AR process described in the question
        f(n)=u(n)-w*u(n-1); % calculate the error in estimation
        w=w+mu*u(n-1)*f(n); % update the weights in the manner associated
                            % LMS based on the recient error in prediction
                            % of the last input
    end
    g=g+f.^2;               % accumulate squared error of estimation
end
g=g/numberOfRuns;           % normalize the accumulated error to reflect
                            % the average error based on the number of
                            % monte Carlo simulations completed
semilogy([1:numberOfDatapoints],g,'k') % plot the MSE vs iterations
title('Graph of squared error vs iteration for Problem 6.16')
xlabel('Number of iterations') % x-axis label
ylabel('Squared Error') % y-axis label