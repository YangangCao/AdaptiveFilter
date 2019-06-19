clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adaptive Filter Theory 5e Solution Manual                          %
%                                                                    %
% Chapter 6                                                          %
% Question 18                                                        %
% Parts a), b)                                                       %
%                                                                    %
% Program written to run on MATLAB 2010a (R)                         %
%                                                                    %
% By Kelvin Hall                                                     %
% June 30, 2014                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
mu=0.01;
H1=[0.25,1,0.25];   % channel characteristics associated with part i)
H2=[0.25,1,-0.25];  % channel characteristics associated with part ii)
H3=[-0.25,1,0.25];  % channel characteristics associated with part iii)
standardDeviationNoise=sqrt(0.01);
delay=22;          % max delay

numberOfDatapoints=400;  % number of iterations per monte carlo run
NumberOfRuns=200;        % number of monte carlo simulations per delay

u=zeros(numberOfDatapoints+delay+1+21,1);   % allocate memory for input to
                                            % channel
                                                   
r=ones(1,numberOfDatapoints+delay+1+21);  % allocate memory for output of
                                          % channel/input to adaptivefilter
                                          
e=ones(numberOfDatapoints+delay+1+21,1); % allocate memory for instantenous
                                         % error
                                         
E=ones(numberOfDatapoints+delay+1+21,22); % allocate memory for cummulative
                                            %mean squared error
for delay=1:22      % index through various delays
    for k=1:NumberOfRuns    % index through monte carlo simulations
        W=zeros(21,1); % rest the weights to zero after each montecarlo run
        
        for n=delay+1+21:numberOfDatapoints+delay+1+21 % iteration
            
            u(n)=rand(1);   % produce a binomial distribution with equal 
            if u(n)>0.5     % probability of recieving a +1 or -1
                u(n)=1;
            else
                u(n)=-1;
            end
            d=u(n-delay);
            r(n)=H1*[u(n);u(n-1);u(n-2)]+randn*standardDeviationNoise;
            e(n)=u(n-delay)-r(n:-1:n-20)*W;
            W=W+mu*e(n)*r(n:-1:n-20)';
        end
        if delay==22
            delay;
        end
        E(:,delay)=E(:,delay)+e.^2;
    end
end
E(:,:)=E(:,:)/NumberOfRuns;
for i=1:22
    for n=1:(numberOfDatapoints+44)
        if E(n,i)==1/NumberOfRuns||E(n,i)==(NumberOfRuns+1)/NumberOfRuns
            E(n,i)=NaN;
        end
    end
end
mesh(E)
plot(E(:,11))
    