%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adaptive Filter Theory 5e Solution Manual                           %
%                                                                     %
% Chapter 6                                                           %
% Question 17                                                         %
% Part a)                                                             %
%                                                                     %
% Program written to run on MATLAB 2010a (R)                          %
%                                                                     %
% By Kelvin Hall                                                      %
% June 30, 2014                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear
LargeNumber=10000;  % as the sample variance approaches the actual variance
                    % of the system as the number of samples grow, the
                    % large nuber is the number of samples taken to
                    % calculate the sample variance.
NoiseVariance=10;    % An initial guess at the nois
NoiseStandardDeviation=sqrt(NoiseVariance);
                    
u=zeros(LargeNumber+3,1); % Allocate memory for input data stream

aOne=0.1;       % AR paramters as given in the textbook
aTwo=-0.8;

Tolerance=0.001;   % the magnitude of the maximal error that is acceptable
MaxNumberOfImpovementCycles=100;  % The maximum number of updating cycles
error=Tolerance+1;  % an initial error larger then the tolerance of so at
                    % least one cycle is completed 
                    
mu=0.02;        % update learning parameter for LMS like learning

k=1;        % number of cycles completed itialized to the value one
while(abs(error)>Tolerance && k<MaxNumberOfImpovementCycles)
    NoiseStandardDeviation=sqrt(NoiseVariance);
    
    for n=3:(LargeNumber+3)
        u(n)=-aOne*u(n-1)-aTwo*u(n-2)+randn(1)*NoiseStandardDeviation;
    end
    
    error=(1-var(u(3:LargeNumber+3))); % Calculate the difference between
                     % process sample variance and desired process variance
                     
    NoiseVariance=NoiseVariance+mu*error; % adjust the noise variance to
                                         % reduce error in process variance
    k=k+1;  % increment loop variable
end

if abs(error)>Tolerance % provide message explaining why solution wasn't found
    fprintf('sorry you may need to provide a closer starting value or')
    fprintf('run the system for more cycles')
    systemVariance=var(u(3:LargeNumber+3))
    NoiseVariance
else   
    SampleProcessVariance=var(u(3:LargeNumber+3))
    NoiseVariance
    k
end