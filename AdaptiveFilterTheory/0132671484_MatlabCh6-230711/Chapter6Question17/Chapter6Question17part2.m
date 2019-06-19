%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adaptive Filter Theory 5e Solution Manual                           %
%                                                                     %
% Chapter 6                                                           %
% Question 17                                                         %
% Parts c), d) and e)                                                 %
%                                                                     %
% Program written to run on MATLAB 2010a (R)                          %
%                                                                     %
% By Kelvin Hall                                                      %
% June 30, 2014                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc
numberOfDatapoints=100; % try 300 to get a better picture of the process
numberOfRuns=100;   % The number of monteCarlo runs necissary for part d)
aOne=0.1;          %AR paramter
aTwo=-0.8;         %AR paramter as per textbook description

NoiseVariance=0.28;% The Noise variance found through part A of the problem
NoiseStandardDeviation=sqrt(NoiseVariance);
mu=0.05;

stream = RandStream('mt19937ar','Seed',30);  % seed the random number
RandStream.setDefaultStream(stream);        % generator for reproducable
                                            % results

u=zeros(numberOfDatapoints+3,1);    % Allocate memory for input data stream
f=zeros(numberOfDatapoints+3,1);    % Allocate memory for error between 
                                    % prediction and resutls
epsilonOne=zeros(numberOfDatapoints+3,1);% Allocate memory for error  
                                   % between found weight and AR paramter 1
epsilonTwo=zeros(numberOfDatapoints+3,1);% Allocate memory for error  
                                   % between found weight and AR paramter 2
g=zeros(numberOfDatapoints+3,1);% allocate memory for squared-averaged error
J=zeros(numberOfDatapoints+3,1);% allocate memory for theoretical error

weights=zeros(2,numberOfDatapoints+3); % allocate memory for weights

for k=1:numberOfRuns % increment through the monte carlo runs of the problem
   W=zeros(2,numberOfDatapoints+3);    % reset the weights to zero between 
                                       % each montecarlo run
   for n=3:numberOfDatapoints+3
        u(n)=aOne*u(n-1)+aTwo*u(n-2)+randn(1)*NoiseStandardDeviation;
             % create a new piece of input data to the filter
             
        f(n)=u(n)-W(1,n-1)*u(n-1)-W(2,n-1)*u(n-2);
             % Caluclate error between desired result and predicted results
             
        epsilonOne(n)=aOne-W(1,n-1); 
             % calculate error between weight one and AR paramter one
             
        epsilonTwo(n)=aTwo-W(2,n-1); 
             % calculate error between weight two and AR paramter two
        
        W(:,n)=W(:,n-1)+mu*f(n)*[u(n-1);u(n-2)]; 
                    %update the weights as per LMS
    end
    g=g+f.^2;      % accumulate squared error of estimation
end
g=g/numberOfRuns;  % accumulate squared error of estimation

for n=3:numberOfDatapoints+3
    J(n)=(1-NoiseVariance*(1+mu/2))*(1-mu)^(2*n)+NoiseVariance*(1+mu/2);
            %create the data for the graph of theoretical values for part
            %e) of the the experiment.
end

index=1:numberOfDatapoints+3; %the x axis of the plot


subplot(2,2,1) % first subplot showing the PSD of the error f, where it can
               % be seen that the error is similar to that of white noise
               % with approximatly even values across the frequency range
plot(index,(abs((fft(f)))).^2)
xlim([0 numberOfDatapoints+3])
title('Power Spectral Plot of f')
xlabel('Iterations') % x-axis label
ylabel('Power Level') % y-axis label

subplot(2,2,2) % second subplot showing the power spectral density of Error
               % epsilon 1, as there is a significant offset error the
               % frequency representation is highly noneven, with more
               % being towards the lower frequency
plot(index,(abs((fft(epsilonOne)))).^2)
xlim([0 numberOfDatapoints+3])
title('Power Spectral Plot of \epsilon_1')
xlabel('Frequency') % x-axis label
ylabel('Power Level') % y-axis label

subplot(2,2,3) % third subplot showing the power spectral density of Error
               % epsilon 2, as there is a significant offset error the
               % frequency representation is highly noneven, with more
               % being towards the lower frequency
plot(index,(abs((fft(epsilonTwo)))).^2)
xlim([0 numberOfDatapoints+3])
title('Power Spectral Plot of \epsilon_2')
xlabel('Iterations') % x-axis label
ylabel('Power Level') % y-axis label

subplot(2,2,4) % Solution to part d) and e) on the same graph for
               % comparison purposes.
plot(index,g(1:numberOfDatapoints+3),index,J)
xlim([0 numberOfDatapoints+3])
ylim([0 1])
legend('Experimental Results','Theoretical Results')
xlabel('Number of iterations') % x-axis label
ylabel('Squared Error') % y-axis label
title({'Comparison Plot of Theoretical Learning Curve';...
                       'to Observed MonteCarlo Results'})
