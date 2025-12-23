% Jonathan Plavnik, CMOR220, Fall 2025, Project 2
% project3_population_model.m
% This project is split into two tasks. The first deals with the logistic
% growth model for the rabbit population in the absence of predactors.
% Using given parameters values, Eulers method is used to approx answers.
% For the second part, we solve the predator prey problem with given
% paramaters and compare eulers method to matlabs built in ode45 solver. We
% compare the populations of prey and predator against time then against
% each other, recording how it changes over time and producing a video.
% Last modified: September 30, 2025


disp("Project 3: Population Model");

% Question Ra

% given parameters
K = 100;
r = 1;
delta = 0.01;
t = 0:delta:15; 

f_derivative = @(R) r*R*(1 - R/K);
%input: evalute the value of the derivative at R

R = zeros(1,length(t)); 
% initial population, given
R(1) = 20; 
for i = 1:length(t)-1
    R(i+1) = R(i) + delta*f_derivative(R(i)); % euler step
end

% plot everything
figure(1);
plot(t,R,'b');
xlabel('Time'); 
ylabel('Rabbit Population');
title('Logistic Population Growth'); 


% Question Rb

% given parameters
K = 100;
r = 1;
delta = 0.01;
t = 0:delta:15; 

f_derivative = @(R) r*R*(1 - R/K); % anonymous func for derivative

R = zeros(1,length(t)); 
% initial population, given. This is the only difference from the code for
% Ra.
R(1) = 150; 
for i = 1:n-1
    R(i+1) = R(i) + delta*f_derivative(R(i)); % calculate next x value (R(i+1)) by multiplying delta x (delta) by delta y (derivative = R(i))
end

figure(2);
plot(t,R,'b');
xlabel('Time'); 
ylabel('Rabbit Population');
title('Logistic Population Growth'); 



% Question Rfa

% initial given parameters
k1 = 3; 
k2 = 3e-3; 
k3 = 6e-4; 
k4 = 0.5;

R0 = 1000; 
F0 = 500;

% create anonymous functions for the derivatives
dR = @(R,F) k1*R - k2*R.*F;
dF = @(R,F) k3*R.*F - k4*F;

% our given delta is .01, we create the basic arrays that we will be using
% to calculate the next y value

% initalize necessary code
delta1 = 0.01;
t = 0:delta1:15; 
R_1 = zeros(1,length(t)); 
F_1 = R_1;
R_1(1)=R0;
F_1(1)=F0;

% euler step
for i=1:length(t)-1
    R_1(i+1) = R_1(i) + delta1*dR(R_1(i),F_1(i));
    F_1(i+1) = F_1(i) + delta1*dF(R_1(i),F_1(i));
end

% same as above but with smaller step of delta .001. Expect slightly more
% accurate results.

delta2 = 0.001;
t_2 = 0:delta2:15; 
R_2 = zeros(1,length(t_2)); 
F_2 = R_2;
R_2(1)=R0;
F_2(1)=F0;
for i=1:length(t)-1
    R_2(i+1) = R_2(i) + delta2*dR(R_2(i),F_2(i));
    F_2(i+1) = F_2(i) + delta2*dF(R_2(i),F_2(i));
end

% RFb

% re state the two differential equations for fox and rabbit populations
predprey = @(t,y)[k1*y(1) - k2*y(1)*y(2);
                  k3*y(1)*y(2) - k4*y(2)];
optionss = odeset('RelTol',1e-6); % set tol
[t_45,Y_45] = ode45(predprey,[0 15],[R0 F0], optionss); % use ode45 to solve


% Rfc

% plot the results from ode45 on one figure to compare results
figure(3);
hold on;
plot(t,R_1,'b--',t,F_1,'r--');
plot(t_2,R_2,'b:',t_2,F_2,'r:');
plot(t_45,Y_45(:,1),'b');
plot(t_45,Y_45(:,2),'r');
xlabel('Time'); 
ylabel('Population');
legend('Rabbits (Euler|=0.01)','Foxes (Euler=0.01)', ...
       'Rabbits (Euler=0.001)','Foxes (Euler=0.001)', ...
       'Rabbits (ODE)','Foxes (ODE)');
title('Predator Prey Population Over Time'); 
hold off;


% answer to questions posed in rfc
% I observed that the aprox sols from eulers method were close to sols from
% ode45. Sol from eulers method using delta = 0.001 was closer than 0.01. The overall
% patterns from the graph tell us that a rise in rabbit population led to
% rise in fox population. When number of foxes is too high rabbits begin to
% die off, and less rabbits leads to less foxes. Less foxes leads to more
% rabbits. This cycle repeats.

% RFd + RFe 

video = VideoWriter('LimitCircle.avi');
open(video);

% necessary steps for showing the video


disp(Y_45);

figure(4);

% plot all the points comparing rabbit pop and fox pop along Y_45 for time
% from t = 0 to t = 15.
for i=1:length(t_45)
    plot(Y_45(1:i,1),Y_45(1:i,2),'-', 'Color', 'w');
    xlabel('Rabbit Population'); 
    ylabel('Fox Population');
    title("Rabbits vs Foxes (Time: " + 15 * i/length(t_45)+ ")"); % change title for time
    axis([0 max(Y_45(:,1)) 0 max(Y_45(:,2))]); % give bounds for the axixes of the graph
    writeVideo(video, getframe(gcf));
end
close(video);

% answer to questions for rfd
% I observd the same relationship that we saw previously when comparing both rabbit
% and foxes populations. The cycle of increasing rabbit population leading
% to increased fox, then rabbit population decreases, then foxes, then
% rabbit increases once again. This relates to the figure we produced in
% part c because instead of plotting both of the populations against time
% we plot against each other. This led to the same pattern being observed.