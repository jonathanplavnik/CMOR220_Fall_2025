% Jonathan Plavnik, CMOR220, Fall 2025, Project 5
% project5_bridge.m
% In this project, we will examine a simple truss bridge and its ability to support the weight of cars
% driving over it. A bridge is defined by its number of sections, having 5
% * nos + 5 fibers and 2 * nos + 2 nodes. We initialize a bridge using a
% function, which gives us helpful matrices to keep track of x, y coords
% and length of each fiber. We then deform this bridge as we simulate the
% affect that cars of varying weights would have if they crossed the
% bridge. We plot all the results at the end and answer specific questions
% about the ability of the bridge.
% Last modified: October 22, 2025


disp("Project 5: Bridge");

% build the adjacency matrix, x and y coords for each of the starting and ending nodes for each of the fibers and length of all the fibers for the basic truss bridge
function [adj, xc, yc, len] = build_basic_bridge(nos)
% list of inputs: nos number of sections
% list of outputs: It returns adj (the adjacency matrix for the bridge), xc (the x-coordinates of the fibers in the bridge), yc (the y-coordinates), and len
% (the length of each fiber).

% calculate some helpful numbers
num_nodes = 2 * nos + 2;
num_fibers = 5 * nos + 5;
s = 1/sqrt(2);
% initialize the return values
adj = zeros(num_fibers, 2*num_nodes);
xc = zeros(num_fibers, 2);
yc = zeros(num_fibers, 2);
len = ones(num_fibers, 1);
% build the left side of bridge
adj(1, 1) = 1;
adj(2, [3 4]) = [s s];
xc(1, :) = [0 1];
xc(2, :) = [0 1];
yc(1, :) = [0 0];
yc(2, :) = [0 1];
len(2) = 1/s;
% build the middle of bridge
for i = 1:nos
    adj(5 * i - 2, [4*i-2 4*i]) = [-1 1];
    adj(5 * i - 1, [4*i-1 4*i 4*i+1 4*i+2]) = [-s s s -s];
    adj(5 * i, [4*i-1 4*i+3]) = [-1 1];
    adj(5 * i + 1, [4*i-3 4*i-2 4*i+3 4*i+4]) = [-s -s s s];
    adj(5 * i + 2, [4*i-3 4*i+1]) = [-1 1];

    xc(5 * i - 2, :) = [i i];
    xc(5 * i - 1, :) = [i i+1];
    xc(5 * i, :) = [i i+1];
    xc(5 * i + 1, :) = [i i+1];
    xc(5 * i + 2, :) = [i i+1];

    yc(5 * i - 2, :) = [0 1];
    yc(5 * i - 1, :) = [1 0];
    yc(5 * i, :) = [1 1];
    yc(5 * i + 1, :) = [0 1];
    yc(5 * i + 2, :) = [0 0];

    len(5 * i - 2, :) = 1;
    len(5 * i - 1, :) = 1/s;
    len(5 * i, :) = 1;
    len(5 * i + 1, :) = 1/s;
    len(5 * i + 2, :) = 1;
end
% build the right side of bridge

adj(5 * nos + 3, [4 * nos + 2 4 * nos + 4]) = [-1 1];
adj(5 * nos + 4, [4 * nos + 3 4 * nos + 4]) = [-s s];
adj(5 * nos + 5, 4 * nos + 1) = -1;

xc(5 * nos + 3, :) = [nos+1 nos+1];
xc(5 * nos + 4, :) = [nos+1 nos+2];
xc(5 * nos + 5, :) = [nos+1 nos+2];

yc(5 * nos + 3, :) = [0 1];
yc(5 * nos + 4, :) = [1 0];
yc(5 * nos + 5, :) = [0 0];

len(5 * nos + 3, :) = 1;
len(5 * nos + 4, :) = 1/s;
len(5 * nos + 5, :) = 1;



end



% deforms the bridge we created using the previous function
function [dx,dy,work,X,Y] = deform_basic_bridge(nos,adj,xc,yc,len,force)
% list of inputs: nos number of sections, adj (the adjacency matrix for the bridge), xc (the x-coordinates of the fibers in the bridge), yc (the y-coordinates), and len
% (the length of each fiber), and force (the column vector associated with the weight being applied to each node)
% list of outputs: dx deformed x coords, dy deformed y coords, work scalar that measures the way the bridge system as a whole has deformed and how much work the
% fibers are doing to prevent this deformation, X and Y are helper matrices that associate each node (not degree of freedom) in the bridge with its displacement in the horizontal and vertical
% directions

% calculate some helpful numbers


stiffness = adj'*diag(1./len)*adj; % This is matrix S. Now we solve Sx = f
displacements = stiffness\force; % This is x
work = displacements'*force; % calculate work
X = displacements(1:2:end); % displacements to  x values of nodes (length is num nodes)
Y = displacements(2:2:end); % displacements to y values of nodes


% initialize the return values
dx = zeros(size(xc));
dy = zeros(size(yc));
% deform the left side of bridge
dx(1,:) = xc(1,:) + [0 X(1)];
dx(2,:) = xc(2,:) + [0 X(2)];
dy(1,:) = yc(1,:) + [0 Y(1)];
dy(2,:) = yc(2,:) + [0 Y(2)];
% deform the middle of bridge
for i = 1:nos
    dx(5 * i - 2, :) = xc(5 * i - 2, :) + [X(2*i-1) X(2*i)];
    dx(5 * i - 1, :) = xc(5 * i - 1, :) + [X(2*i) X(2*i+1)];
    dx(5 * i, :) = xc(5 * i, :) + [X(2*i) X(2*i+2)];
    dx(5 * i + 1, :) = xc(5 * i + 1, :) + [X(2*i-1) X(2*i+2)];
    dx(5 * i + 2, :) = xc(5 * i + 2, :) + [X(2*i-1) X(2*i+1)];

    dy(5 * i - 2, :) = yc(5 * i - 2, :) + [Y(2*i-1) Y(2*i)];
    dy(5 * i - 1, :) = yc(5 * i - 1, :) + [Y(2*i) Y(2*i+1)];
    dy(5 * i, :) = yc(5 * i, :) + [Y(2*i) Y(2*i+2)];
    dy(5 * i + 1, :) = yc(5 * i + 1, :) + [Y(2*i-1) Y(2*i+2)];
    dy(5 * i + 2, :) = yc(5 * i + 2, :) + [Y(2*i-1) Y(2*i+1)];
    
end
% deform the right side of bridge
dx(end - 2, :) = xc(end-2, :) + [X(2*nos + 1) X(2*nos + 2)];
dx(end - 1, :) = xc(end-1, :) + [X(2*nos + 2) 0];
dx(end, :) = xc(end, :) + [X(2*nos + 1) 0];

dy(end - 2, :) = yc(end-2, :) + [Y(2*nos + 1) Y(2*nos + 2)];
dy(end - 1, :) = yc(end-1, :) + [Y(2*nos + 2) 0];
dy(end, :) = yc(end, :) + [Y(2*nos + 1) 0];

% return
end

% build a basic bridge, then load it with the weight of several cars, deform the bridge, 
% and then call plot bridge to plot the results.

function build_load_plot_basic_bridge(nos, car_weight)
% list of inputs: nos number of sections, car_weight weight of the car
% which will give us the force vector
% list of outputs: No return outputs but will create several plots

% initiailize bridge
[adj, xc, yc, len] = build_basic_bridge(nos);
if(car_weight == 0)
    % plot spy
    plot_spy(adj, nos);
    % plot basic bridge no weight thus no deformation
    figure;
    plot_bridge(nos, xc, yc);
    title("Bridge with " + nos + " section(s)", 'FontWeight','bold');
    h = findobj(gca,'Type','line');
    set(h,'LineWidth',3,'Color',[0.30 0.50 1]);
else
    % initialize force
    force = zeros(2*(2*nos + 2), 1);
    force(2:4:end) = -car_weight;
    
    [dx,dy,work,X,Y] = deform_basic_bridge(nos,adj,xc,yc,len,force);
    
    % plot the bridge
    figure;
    plot_bridge(nos, dx, dy);
    title(nos + " Section Deformed Bridge", 'FontWeight','bold');
    subtitle("Car Weights = " + car_weight + " | Work = " + work);
    h = findobj(gca,'Type','line');
    set(h, 'LineWidth', 3, 'Color', [0.75 0.10 0.15]);
end

end

% create pretty plots of bridges given the inputs
function plot_bridge(nos, xc, yc)
% list of inputs: nos number of sections, xc x coordinates of starting and
% ending nodes for each fiber, yc y coords of starting and ending nodes for
% each fiber
% list of outputs: no return outputs but will plot bridge

% transpose to graph
xc = xc';
yc = yc';
hold on;
axis equal
% change color of axes to not be black since that will mess with fill
ax = gca;
ax.Color = [0.95 0.95 0.95];

xlim([-1, nos + 3]); % set the limits of the graph for clarity
ylim([-1, 1.1]);


fill([-1, 0, .5, -1],[0, 0, -1, -1],'k') % left side black fill
fill([2+nos, 3+nos, 3+nos, 1.5+nos],[0, 0, -1, -1],'k') % right side black fill

%plot the x and y coords
plot(xc, yc, '-', 'LineWidth', 3, 'Color', 'r');
hold off;
end

% create spy matrix for adj matrix
function plot_spy(adj, nos)
% list of inputs: adj adjacency matrix, nos number of sections
% list of outputs: returns nothing but plots spy for adj matrix
figure;
spy(adj);
title("Adjacency Matrix for " + nos + "-Section Bridge", 'FontWeight','bold');
xlabel("Nonzero Matrix Columns");
ylabel("Nonzero Matrix Rows");
end


% main

%plot all 12 requiered graphs
build_load_plot_basic_bridge(1, 0);
build_load_plot_basic_bridge(1, 0.01);
build_load_plot_basic_bridge(1, 0.05);
build_load_plot_basic_bridge(2, 0);
build_load_plot_basic_bridge(2, 0.01);
build_load_plot_basic_bridge(2, 0.05);
build_load_plot_basic_bridge(3, 0);
build_load_plot_basic_bridge(3, 0.01);
build_load_plot_basic_bridge(3, 0.05);



% Response to questions

% 1) If you are driving a light-weight car (0.01 units of weight), what is the maximum length (in
% nos) the bridge could be while still being considered safe? (What is your criterion for “safe”?)
%  First we must define our tolerance for deformation to any node of the bridge. This can be decided somewhat arbitrarily, 
% but I will choose 0.2 units of deformation based on visual analysis of
% the plots for 0.01 units of weight. Using this definition, I examined
% what happened to the bridge under differing nos. I gathered that while
% the absolute value of all elements of X and Y is under 0.2 for nos = 1,
% 2,3, this is not the case for nos = 4. Therefore I have determined the
% maximum nos to be considered safe is 3. Consider the code below to prove
% this outcome.

for(nos = 3:4)
    [adj, xc, yc, len] = build_basic_bridge(nos);
    force = zeros(2*(2*nos + 2), 1);
    force(2:4:end) = -.01;
    [dx,dy,work,X,Y] = deform_basic_bridge(nos,adj,xc,yc,len,force);
    disp("Nos = " + nos + ", Value of X (displacements");
    disp(X);
    disp("Nos = " + nos + ", Value of Y (displacements");
    disp(Y);
end


% 2. What about if you are driving an 18-wheeler truck (0.05 units of weight)?
% Applying the same definition as last time, we once again consider the
% deformation values that can be found in X, Y. We are looking to ensure
% that none of the entries in these vectors is greater than 0.2. The only
% value of nos for which this does not occur is nos = 1, therefore only one
% section is safe Consider the code below to prove
% this outcome.

for(nos = 1:3)
    [adj, xc, yc, len] = build_basic_bridge(nos);
    force = zeros(2*(2*nos + 2), 1);
    force(2:4:end) = -.01;
    [dx,dy,work,X,Y] = deform_basic_bridge(nos,adj,xc,yc,len,force);
    disp("Nos = " + nos + ", Value of X (displacements");
    disp(X);
    disp("Nos = " + nos + ", Value of Y (displacements");
    disp(Y);
end

% 3. How does the spy of the adjacency matrix A for a particular bridge relate to how it deforms?
% The spy plot of adj matrix A represents the number of nonzero elements,
% indicating some kind of connection between nodes. We have termed this to
% be a fiber. These fibers will change under deformation hence these vary, or in other words 
% they represent the degrees of freedom of the bridge. The greater these vary the greater deformation occurs,
% thus bridges with higher degrees of freedom will deform more.

% 4. What trends do you notice as nos increases? (Relate your answer to the spy of a particular
%bridge, how the bridge deforms under different stress conditions, etc.)

% As the number of sections increase, I noticed the number of fibers represented in 
% the spy of adj matrix increasing. As per the correlation identified in
% question 3, this will lead to more deformation. Longer bridges will
% deform more. To better illustrate this point, I have repeated the plots
% of bridges under 0.05 units of weight, contrasting bridges with 1 and 3
% sections. They are written below.
%build_load_plot_basic_bridge(1, 0.05);
%build_load_plot_basic_bridge(3, 0.05);
