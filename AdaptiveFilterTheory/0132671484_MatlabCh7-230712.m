%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adaptive Filter Theory 5e Solution Manual                              %
%                                                                        %
% Chapter 7                                                              %
% Question 10                                                            %
% Parts a), b), and c)                                                   %
%                                                                        %
% Program written to run on MATLAB 2010a (R)                             %
%                                                                        %
% By Kelvin Hall                                                         %
% July 2, 2014                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc
numberOfDatapoints=100; % try 300 to get a better picture of the process
numberOfRuns=100;   % The number of monteCarlo runs necissary for part d)
aOne=0.1;          %AR paramter
aTwo=-0.8;         %AR paramter as per textbook description

NoiseVariance=0.28;% The Noise variance found through part A of the problem
NoiseStandardDeviation=sqrt(NoiseVariance);
mu=0.2;
delta=[0.5,0.25,0.75];

u=zeros(numberOfDatapoints+3,1);    % Allocate memory for input data stream
f=zeros(numberOfDatapoints+3,1);    % Allocate memory for error between 
                                    % prediction and resutls
epsilonOne=zeros(numberOfDatapoints+3,1);% Allocate memory for error  
                                   % between found weight and AR paramter 1
epsilonTwo=zeros(numberOfDatapoints+3,1);% Allocate memory for error  
                                   % between found weight and AR paramter 2
J=zeros(numberOfDatapoints+3,1);% allocate memory for theoretical error

weights=zeros(2,numberOfDatapoints+3); % allocate memory for weights

for ell=1:3
    
stream = RandStream('mt19937ar','Seed',2);  % seed the random number
RandStream.setDefaultStream(stream);        % generator for reproducable
                                            % results
g=zeros(numberOfDatapoints+3,1);% allocate memory for squared-averaged error
for k=1:numberOfRuns % increment through the monte carlo runs of the problem
   W=zeros(2,numberOfDatapoints+3);    % reset the weights to zero between 
                                       % each montecarlo run
   for n=3:numberOfDatapoints+3
        u(n)=aOne*u(n-1)+aTwo*u(n-2)+randn(1)*NoiseStandardDeviation;
             % create a new piece of input data to the filter
        
        f(n)=u(n)-W(1,n-1)*u(n-1)-W(2,n-1)*u(n-2);
             % Caluclate error between desired result and predicted results
        
        W(:,n)=W(:,n-1)+(mu/(delta(ell)+u(n-1)^2+u(n-2)^2)*([u(n-1);u(n-2)]*f(n))); 
                    %update the weights as per LMS
        e(:,n)=W(:,n)-[u(n-1);u(n-2)];
    end
    g=g+f.^2;      % accumulate squared error of estimation
end
g=g/numberOfRuns;  % accumulate squared error of estimation
index=1:numberOfDatapoints+3; %the x axis of the plot

subplot(3,3,1+(ell-1)*3)
plot(index,W(1,:))
xlim([0 numberOfDatapoints+3])
ylim([-0.1, 0.5])
str=sprintf('error in tap weight 1 estimate for delta = %.2f',delta(ell));
title(str)
xlabel('Iterations')
ylabel('error between weight and AR paramter') 

subplot(3,3,2+(ell-1)*3) 
plot(index,W(2,:))
xlim([0 numberOfDatapoints+3])
str=sprintf('error in tap weight 2 estimate for delta = %.2f',delta(ell));
title(str)
xlabel('Iterations') 
ylabel('error between weight and AR paramter') 

subplot(3,3,3+(ell-1)*3) 
plot(index,g(1:numberOfDatapoints+3))
xlim([0 numberOfDatapoints+3])
ylim([0 1])
xlabel('Number of iterations')
ylabel('Squared Error')
str=sprintf('Plot of Normalized LMS learning curve with delta =%.2f',delta(ell));
title(str);
end