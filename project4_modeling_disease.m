% Jonathan Plavnik, CMOR220, Fall 2025, Project 4
% project4_modeling_diseases.m
% This project models the SIR model given 3 different set of assumptions.
% Initially part 1 considers the simpleSIR model where total population
% remains constant, part 2 considers the variableSIR model where total
% population varies while deaths not from infection and balanced by births
% and part 3 considers the the model from part 2 but there is a loss of
% immunity meaning that recovered have a fixed percent chance of getting
% sick again. In each of the parts the population of each of these groups
% is tracked day by day and plotted over time and in parts 2 and 3, the
% infected population is modeled as a function of the susceptible function.
% Last modified: October 9, 2025 (Extension)


disp("Project 4: Modeling Infectious Diseases");


% Part 1: SIR with constant total population

% S_n is the number of susceptible individuals at day n, (healthy)
% I_n is the number of infected individuals at day n, (have virus, infect others) 
% R_n is the number of recovered individuals at day n, (had virus, cannot get it, cannot infect)

% Part 1 assumption: S_n + I_n + R_n = M for all n >= 0

% assume that each infected individual has a fixed number α of contacts per day that are sufficient to spread the disease
% assume that a fixed fraction β of the infected group will recover during any given day


% this function will return vectors containing S_n, I_n, R_n for the simple
% SIR model.
function [Sval, Ival, Rval] = simpleSIR(M, alpha, beta, initialval, Tfinal)
% list of inputs: total population M, each infected individual has a fixed
% number α of contacts per day that are sufficient to spread the disease, a
% fixed fraction β of the infected group will recover during any given day,
% the initial values for Sval, Rval, and the number of days the model will
% run for.
% list of outputs: The number of susceptbile, infected, and recovered
% individuals at days 1 through n.



% create empty vectors that will serve as output vectors
Sval = zeros(1, Tfinal + 1);
Ival = Sval;
Rval = Sval;
t = 0:Tfinal;
% initialization of the variables
Sval(1) = initialval(1);
Rval(1) = initialval(2);
Ival(1) = M-Sval(1)-Rval(1);
% loop over the time steps
for i=1:Tfinal
    Sval(i+1)=Sval(i)-((alpha/M)*Sval(i)*Ival(i));
    Rval(i+1)=Rval(i)+(beta*Ival(i));
    Ival(i+1) = M-Sval(i+1)-Rval(i+1);
end

% plot the evolution of the three population groups over 150 days
figure(1);
hold on;
plot(t, Sval, 'b');
plot(t, Rval, 'r');
plot(t, Ival, 'y');
xlabel('nb of days'); 
ylabel('population');
title('Figure 1 (SIR Model Population Constant over 150 days)');
legend('Susceptibles','Recovered','Infectious');
hold off;

end

% calling the function
S_0 = 7.9e6 - 10;
R_0 = 0;
[S, I, R] = simpleSIR(7.9e6, .7, .1, [S_0, R_0], 150);

% Answer the question
% The reason why Sn, In, and Rn changed so significantly with these inputs
% is because alpha, the rate of infection, is high at 0.7 while beta, the rate of
% recovery, is low at 0.1. This mean that In increases sharply at the
% beginning while Sn crashes. Then the high number of In means that despite
% a lower beta, the number of recovered Rn continues to grow until the
% entire population has gotten sick then recovered.




% Part 2: SIR with variable total stable population


% this function will return vectors containing S_n, I_n, R_n for the
% variable SIR model, meaning total population is not fixed and the number
% of deaths not infection is equal to number of births
function [Sval, Ival, Rval] = variableSIR(alpha, beta, gamma, mu, initialval, Tfinal)
% list of inputs: total population M, each infected individual has a fixed
% number α of contacts per day that are sufficient to spread the disease, a
% fixed fraction β of the infected group will recover during any given day,
% Death due to infection will cause a loss of individuals from the infected
% group at a rate gamma, deaths from all over causes at a rate mu, the 
% initial values for Sval, Rval, and the number of days the model will
% run for.
% list of outputs: The number of susceptbile, infected, and recovered
% individuals at days 1 through n.

% create empty vectors that will serve as output vectors
Sval = zeros(1, Tfinal + 1);
Ival = Sval;
Rval = Sval;
M = Sval;
t = 0:Tfinal;

% initialization of the variables
Sval(1) = initialval(1);
Rval(1) = initialval(2);
Ival(1) = initialval(3);
M(1) = Sval(1) + Ival(1) + Rval(1);

% loop over the time steps

for i=1:Tfinal
    Sval(i+1)=Sval(i)-((alpha/M(i))*Sval(i)*Ival(i)) + mu * M(i) - mu * Sval(i);
    Rval(i+1)=Rval(i)+(beta*Ival(i)) - mu * Rval(i);
    Ival(i+1)=Ival(i)+((alpha/M(i))*Sval(i)*Ival(i))- (beta*Ival(i)) - (mu + gamma) * Ival(i);
    M(i+1) = Sval(i+1) + Ival(i+1) + Rval(i+1);
end

% Plot the figure to answer 2.1 (the evolution of each group and the total
% population over time)
figure(2);
hold on;
plot(t, Sval, 'b');
plot(t, Rval, 'r');
plot(t, Ival, 'y');
plot(t, M, 'g');
xlabel('nb of days'); 
ylabel('population');
title('Figure 2 (SIR Model Population Variable & Stable over 1460 days)');
legend('Susceptibles','Recovered','Infectious', 'Total Population');
hold off;

% Plot the figure to answer 2.2 (the plot of infected as a function of susceptibles)
figure(3);
hold on;
plot(Sval, Ival)
xlabel('population of susceptibles'); 
ylabel('population of infected');
title('Figure 3: Infected Population as a function of Susceptible Population (Based on Figure 2)');
hold off;

end


% call the function
[Sn, In, Rn] = variableSIR(.5, 1/3, .01, 1/(76*365), [7.9e6 - 10, 0, 10], 4*365);



% Part 3: SIR with variable total population and loss of immunity


% new assumption: fraction omega of Rn moves to Sn+1



% this function will return vectors containing S_n, I_n, R_n for the variable and loss of immunity SIR model, meaning total population is not fixed and the number
% of deaths not infection is equal to number of births and people who have
% recovered can get sick again
function [Sval, Ival, Rval] = variableimmSIR(alpha, beta, gamma, mu, omega, initialval, Tfinal)
% list of inputs: total population M, each infected individual has a fixed
% number α of contacts per day that are sufficient to spread the disease, a
% fixed fraction β of the infected group will recover during any given day,
% Death due to infection will cause a loss of individuals from the infected
% group at a rate gamma, deaths from all over causes at a rate mu, A fraction ω of
% the recovered population moves to the group of susceptibles, the 
% initial values for Sval, Rval, and the number of days the model will
% run for.
% list of outputs: The number of susceptbile, infected, and recovered
% individuals at days 1 through n.


% create empty vectors that will serve as output vectors
Sval = zeros(1, Tfinal + 1);
Ival = Sval;
Rval = Sval;
M = Sval;
t = 0:Tfinal;

% initialization of the variables
Sval(1) = initialval(1);
Rval(1) = initialval(2);
Ival(1) = initialval(3);
M(1) = Sval(1) + Ival(1) + Rval(1);

% loop over the time steps

for i=1:Tfinal
    Sval(i+1)=Sval(i)-((alpha/M(i))*Sval(i)*Ival(i)) + mu * M(i) - mu * Sval(i) + omega * Rval(i);
    Rval(i+1)=Rval(i)+(beta*Ival(i)) - mu * Rval(i) - omega * Rval(i);
    Ival(i+1)=Ival(i)+((alpha/M(i))*Sval(i)*Ival(i))- (beta*Ival(i)) - (mu + gamma) * Ival(i);
    M(i+1) = Sval(i+1) + Ival(i+1) + Rval(i+1);
end


% plot to question 3.1 (the evolution of each group and the total
% population over time)
figure(4);
hold on;
plot(t, Sval, 'b');
plot(t, Rval, 'r');
plot(t, Ival, 'y');
plot(t, M, 'g');
xlabel('nb of days'); 
ylabel('population');
title('Figure 4 (SIR Model Population Variable, Stable, Loss of Immunity over 1460 days)');
legend('Susceptibles','Recovered','Infectious', 'Total Population');
hold off;


% Plot the figure to answer 3.2 (the plot of infected as a function of susceptibles)

figure(5);
hold on;
plot(Sval, Ival)
xlabel('population of susceptibles'); 
ylabel('population of infected');
title('Figure 5: Infected Population as a function of Susceptible Population (Based on Figure 4)');
hold off;

end


% call the function with given inputs
[Sn, In, Rn] = variableimmSIR(.5, 1/3, .01, 1/(76*365), 1/365, [7.9e6 - 10, 0, 10], 4*365);

