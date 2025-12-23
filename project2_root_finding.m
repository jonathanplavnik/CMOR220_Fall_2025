% Jonathan Plavnik, CMOR220, Fall 2025, Project 2
% project2_root_funding.m
% Find the decay rate (the rate at which the heat dissipates at the end of the bar is the square of the least,
% strictly positive solution, x) using the bisection method and plot the
% decay rate against the length of the bar. In the second part, compare the
% speed of bisection and newton methods to find the root for tolerances
% ranging from 10^-8 to 10^-1.
% Last modified: September 23, 2025


disp("Project 2: Root Finding Methods");

% Newton's method , calls coolfun and coolfundx
function [x , iter] = denewt(x , tol ,L)
% list of inputs: initial guess x, tolerance tol, and length of bar L
% list of outputs: final approximation of root x, number of iterations it
% took
iter = 0;
% while loop to iterate over new guesses using taylor series approximation
% to find approx root close enough to meet tol condition
while(abs(coolfun(x,L)) > tol)
    x_new = x - coolfun(x, L)/coolfundx(x, L); % iterative step to find new guess
    x = x_new;
    iter = iter + 1;
end
end

% for a given a , b , t , and L find x: code bisection method
function [x , iter ] = debis(a ,b ,t ,L)
% list of inputs: left side guess a, right side guess b, tolerance t,
% length of bar L
% list of outputs: final approximation of root x, number of iterations it
% took

iter = 0;
x = (a + b )/ 2;
if(coolfun(a, L) * coolfun(b, L) > 0) % check for satisfactory conditions
    disp("Bisection Method Exit: Solution does not exist in " + a + ", " + b);
elseif(coolfun(a, L) * coolfun(b, L) == 0) % check for edge case when a or b are roots
    if(coolfun(a, L) == 0)
        x = a;
    else
        x = b;
    end
else
    while(abs(coolfun(x, L)) > t) % loop to find root approx, implement bisection using mid point as new guess
        if coolfun(a, L) * coolfun(x, L) < 0
            b = x;
        else
            a = x;
        end
        x = (a + b) / 2; % update guess
        iter = iter + 1;
    end
end
end

% for a given x and L evaluate " cool "
function val = coolfun(x ,L)
% list of inputs: variable input x, length of bar L
% list of outputs: evaluate the following function for the input x
val = sin(x*L) + x *cos(x*L);
end

function val = coolfundx(x ,L)
% list of inputs: variable input x, length of bar L
% list of outputs: evaluate the derivative of the function in coolfun given
% input x
val = L * (cos(x * L) + cos(x* L) - L * x * sin(x * L));
% evaluate the derivative , with respect to x , of coolfun
end

% main code
% set a , t and a range of L and b and find x
% and plot x ^2 against L
% calls debis and denewt and plots their iters per tol

% create a vector for L
L = 1 : 0.1 : 4;
for cnt = 1: length(L)
    a = 0.1;
    t = 0.01;
    [x(cnt), ~] = debis(a, 3/L(cnt), t, L(cnt)); % for each length L generate root and save all
end


% plot all the root approximations on one plot
figure;
plot(L , x .^2, '-x');
title("Cooling Rate vs. Bar Length");
xlabel("Length, L");
ylabel("Decay Rate, x^2");


% given data
L = 1;
a = 0.1;
b = 3;
x = (a+b)/2;
t = 10.^(-8:-1);
for i = 1:length(t) % run newton and bisection for each of the tolerances
    [~, bisectionIteration(i)] = debis(a, b, t(i), L);
    [~, newtonIteration(i)] = denewt(x, t(i), L);
end

% plot the comparison on the same graph using semilogx
figure;
semilogx(t, bisectionIteration, '-x', t, newtonIteration, '-o');
hold on
legend(["Bisection", "Newton"]);
title("Bisection vs. Newton");
xlabel("Tolerance");
ylabel("# of Tterations");