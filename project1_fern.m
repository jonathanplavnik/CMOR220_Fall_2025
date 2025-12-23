% Jonathan Plavnik, CMOR220, Fall 2025, Project 1
% project1_fern.m
% Create a stunning fern plot twice by employing an iterative process along
% with matrix multiplication. Depending on probability, starting vector
% z=[0, 1] will be modified by a predetermined amount. There are more cases
% for the advanced fern. This iteration occurs 4000 times before the final
% fern is formed on one figure, preformed once for each of the ferns.
% Last modified: September 14, 2025


disp("Project 1: Grow a Fern")
z = [1;0];
figure;
hold on;
% We have established the starting conditions.


% Now we iterate 4000 times, with a 39.94% probability of z = [.4, -0.3733;
% .06, .6] * z + [.3533, 0];, and a 60.06% probability of z = [-.8,
% -0.1867; .1371, .8] * z + [1.1, .1];, while plotting each point.
for i = 1:4000
    r = rand;
    if r < .3994
        z = [.4, -0.3733; .06, .6] * z + [.3533;0];
    else
        z = [-.8, -0.1867; .1371, .8] * z + [1.1;.1];
    end
    plot(z(1),z(2),'.','MarkerSize',5);
end

% We give title and axis labels.
title("Simple Fern");
xlabel("X axis");
ylabel("Y axis");
hold off;


z = [1;0];
figure;
hold on;
% We have established the starting conditions.

% Now we iterate 4000 times, with certain probabilites affecting the transformation to z while plotting each point on the same figure.
for i = 1:4000
    r = rand;
    if r < 0.01
        z = [0 0; 0 0.16] * z;
    elseif r < 0.76
        z = [0.85 0.04; -0.04 0.85] * z + [0;1.6];
    elseif r < 0.88
        z= [0.2 -0.26; 0.23 0.22] * z + [0;1.6];
    else
        z = [-0.15 0.28; 0.26 0.24] * z + [0;0.44];
    end
    plot(z(1),z(2),'.','MarkerSize',5);
end

% We give title and axis labels.
title("Advanced Fern");
xlabel("X axis");
ylabel("Y axis");
hold off;