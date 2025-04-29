% script_fcn_VSkel_polytopeShrinkFromEdges
% Tests function: fcn_VSkel_polytopeShrinkFromEdges

% REVISION HISTORY:
% 2021_08_02
% -- first written by S. Brennan within MapGen library
% 2025_04_29 by Sean Brennan
% -- updated main script for VSkel library implementation

close all;


%% Demonstration Examples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  _____                                 _             _   _               ______                           _
% |  __ \                               | |           | | (_)             |  ____|                         | |
% | |  | | ___ _ __ ___   ___  _ __  ___| |_ _ __ __ _| |_ _  ___  _ __   | |__  __  ____ _ _ __ ___  _ __ | | ___  ___
% | |  | |/ _ \ '_ ` _ \ / _ \| '_ \/ __| __| '__/ _` | __| |/ _ \| '_ \  |  __| \ \/ / _` | '_ ` _ \| '_ \| |/ _ \/ __|
% | |__| |  __/ | | | | | (_) | | | \__ \ |_| | | (_| | |_| | (_) | | | | | |____ >  < (_| | | | | | | |_) | |  __/\__ \
% |_____/ \___|_| |_| |_|\___/|_| |_|___/\__|_|  \__,_|\__|_|\___/|_| |_| |______/_/\_\__,_|_| |_| |_| .__/|_|\___||___/
%                                                                                                    | |
%                                                                                                    |_|
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Demonstration%20Examples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Demonstration case 1: non-normal wall shrinking
fig_num = 9001;
figure(fig_num);
clf;

% this polytope has a vertical wall
vertices = [0 0; 2 0; 1 2; 0 1; 0 0]*5;
% assert that the wall is vertical to start (i.e. x position of 1st and 4th
% vertices are equal)
assert(vertices(1,1) == vertices(4,1));
test_polytope.vertices = vertices;
% fill in other fields from the vertices field
test_polytope = fcn_VSkel_fillPolytopeFieldsFromVertices(test_polytope);


% perform a small edge shrink
edge_cut = 0.1;
[shrunk_polytope, new_vertices, new_projection_vectors, cut_distance] = fcn_VSkel_polytopeShrinkFromEdges(test_polytope, edge_cut,fig_num);

% Check variable types
assert(isstruct(shrunk_polytope));
assert(iscell(new_vertices));
assert(iscell(new_projection_vectors));
assert(isnumeric(cut_distance));

% Check variable sizes
num_nested = 3; % This is the number of nested figures within the vertex skeleton
assert(length(new_vertices)==num_nested);
assert(length(new_projection_vectors)==num_nested);
assert(isscalar(cut_distance(1,:)));
assert(length(cut_distance(:,1))==num_nested);

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));


% extract new vertices
new_vertices = shrunk_polytope.vertices;
% assert that new vertices are within 5% error of having the same x
% position
error_tolerance = 0.05;
vertical_error = abs(new_vertices(1,1)-new_vertices(4,1))/new_vertices(4,1);
% note this error is currently about 33% as of the commit this line was
% added
assert(vertical_error <= error_tolerance,['Wall should be vertical but x positions of start and end point',...
    ' were %d and %d yielding a vertical error of %d. Error tolerance was %d.\n'],...
    new_vertices(1,1),new_vertices(4,1),vertical_error,error_tolerance);

%% failing case 2: concave polytope produced
fig_num = 9002;
figure(fig_num);
clf;

% this polytope is convex
vertices = [0 0; 2 0; 1 2; 0 1; 0 0]*5;
test_polytope.vertices = vertices;
% fill in other fields from the vertices field
test_polytope = fcn_VSkel_fillPolytopeFieldsFromVertices(test_polytope);
% assert that the polytope is convex to start
[angles, ~, ~] = fcn_VSkel_polytopeFindVertexAngles(test_polytope.vertices);
interior_angles = 180-angles*180/pi;
assert(~any(interior_angles>180));
% perform a large edge shrink
edge_cut = 3.6;

[shrunk_polytope, new_vertices, new_projection_vectors, cut_distance] = fcn_VSkel_polytopeShrinkFromEdges(...
    test_polytope,edge_cut,fig_num);


% Check variable types
assert(isstruct(shrunk_polytope));
assert(iscell(new_vertices));
assert(iscell(new_projection_vectors));
assert(isnumeric(cut_distance));

% Check variable sizes
num_nested = 3; % This is the number of nested figures within the vertex skeleton
assert(length(new_vertices)==num_nested);
assert(length(new_projection_vectors)==num_nested);
assert(isscalar(cut_distance(1,:)));
assert(length(cut_distance(:,1))==num_nested);

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

% extract new vertices
new_vertices = shrunk_polytope.vertices;
% assert that the polytope is convex after shrinking
[new_angles, ~, ~] = fcn_VSkel_polytopeFindVertexAngles(new_vertices);
new_interior_angles = 180-new_angles*180/pi;
assert(~any(new_interior_angles>180),['All interior angles must be < 180 ',...
    'polytope to be convex']);

%% failing case 3: NaN angles produced
fig_num = 9003;
figure(fig_num);
clf;

% this polytope is convex
vertices = [0 0; 2 0; 1 2; 0 1; 0 0]*5;
test_polytope.vertices = vertices;
% fill in other fields from the vertices field
test_polytope = fcn_VSkel_fillPolytopeFieldsFromVertices(test_polytope);
% assert that the polytope is convex to start
[angles, ~, ~] = fcn_VSkel_polytopeFindVertexAngles(test_polytope.vertices);
interior_angles = 180-angles*180/pi;
assert(~any(interior_angles>180));

% perform an even larger edge shrink - this is larger than the object is
edge_cut = 4;
[shrunk_polytope, new_vertices, new_projection_vectors, cut_distance] = fcn_VSkel_polytopeShrinkFromEdges(...
    test_polytope,edge_cut,fig_num);


% Check variable types
assert(isstruct(shrunk_polytope));
assert(iscell(new_vertices));
assert(iscell(new_projection_vectors));
assert(isnumeric(cut_distance));

% Check variable sizes
num_nested = 3; % This is the number of nested figures within the vertex skeleton
assert(length(new_vertices)==num_nested);
assert(length(new_projection_vectors)==num_nested);
assert(isscalar(cut_distance(1,:)));
assert(length(cut_distance(:,1))==num_nested);

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Repeated cuts
fig_num = 9004;
figure(fig_num);
clf;


vertices = [0 0; 0.4 0.1; 1 1; 0 1; 0 0]*5;
test_polytope.vertices = vertices;

% fill in other fields from the vertices field
test_polytope = fcn_VSkel_fillPolytopeFieldsFromVertices(test_polytope);

step =  0.1;
for edge_cut = step:step:2
    [shrunk_polytope, new_vertices, new_projection_vectors, cut_distance] = fcn_VSkel_polytopeShrinkFromEdges(...
        test_polytope,edge_cut,fig_num);
end


% Check variable types
assert(isstruct(shrunk_polytope));
assert(iscell(new_vertices));
assert(iscell(new_projection_vectors));
assert(isnumeric(cut_distance));

% Check variable sizes
num_nested = 3; % This is the number of nested figures within the vertex skeleton
assert(length(new_vertices)==num_nested);
assert(length(new_projection_vectors)==num_nested);
assert(isscalar(cut_distance(1,:)));
assert(length(cut_distance(:,1))==num_nested);

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Random polytope calculation
fig_num = 9005;
figure(fig_num);
clf;

rng(1);

% Set up polytopes
load('testData_fcn_VSkel_plotPolytopes.mat','polytopes','polytopes2');

% Pick a random polytope
Npolys = length(polytopes);
rand_poly = 1+floor(rand*Npolys);
shrinker = polytopes(rand_poly);

% for edge_cut = edge_cut_step:edge_cut_step:(shrinker.max_radius+edge_cut_step)
edge_cut_step = 0.01;
for edge_cut = edge_cut_step:edge_cut_step:(shrinker.max_radius/1.5+edge_cut_step)
    [shrunk_polytope, new_vertices, new_projection_vectors, cut_distance] = fcn_VSkel_polytopeShrinkFromEdges(...
        shrinker,edge_cut,fig_num);
end

% Check variable types
assert(isstruct(shrunk_polytope));
assert(iscell(new_vertices));
assert(iscell(new_projection_vectors));
assert(isnumeric(cut_distance));

% Check variable sizes
num_nested = 5; % This is the number of nested figures within the vertex skeleton
assert(length(new_vertices)==num_nested);
assert(length(new_projection_vectors)==num_nested);
assert(isscalar(cut_distance(1,:)));
assert(length(cut_distance(:,1))==num_nested);

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Basic testing examples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  ____            _        _______        _   _                ______                           _
% |  _ \          (_)      |__   __|      | | (_)              |  ____|                         | |
% | |_) | __ _ ___ _  ___     | | ___  ___| |_ _ _ __   __ _   | |__  __  ____ _ _ __ ___  _ __ | | ___  ___
% |  _ < / _` / __| |/ __|    | |/ _ \/ __| __| | '_ \ / _` |  |  __| \ \/ / _` | '_ ` _ \| '_ \| |/ _ \/ __|
% | |_) | (_| \__ \ | (__     | |  __/\__ \ |_| | | | | (_| |  | |____ >  < (_| | | | | | | |_) | |  __/\__ \
% |____/ \__,_|___/_|\___|    |_|\___||___/\__|_|_| |_|\__, |  |______/_/\_\__,_|_| |_| |_| .__/|_|\___||___/
%                                                       __/ |                             | |
%                                                      |___/                              |_|
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Basic%20Testing%20%20Examples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Basic example of vertex calculation - a square
fig_num = 1001;
figure(fig_num);
clf;

vertices = [0 0; 1 0; 1 1; 0 1; 0 0]*5;
% vertices = [0 0; 2 0; 1 2; 0 1; 0 0]*5;
test_polytope.vertices = vertices;

% fill in other fields from the vertices field
test_polytope = fcn_VSkel_fillPolytopeFieldsFromVertices(test_polytope);


edge_cut = 0.1;
[shrunk_polytope, new_vertices, new_projection_vectors, cut_distance] = fcn_VSkel_polytopeShrinkFromEdges( test_polytope,edge_cut,fig_num);

% Check variable types
assert(isstruct(shrunk_polytope));
assert(iscell(new_vertices));
assert(iscell(new_projection_vectors));
assert(isnumeric(cut_distance));

% Check variable sizes
num_nested = 2; % This is the number of nested figures within the vertex skeleton
assert(length(new_vertices)==num_nested);
assert(length(new_projection_vectors)==num_nested);
assert(isscalar(cut_distance(1,:)));
assert(length(cut_distance(:,1))==num_nested);

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Basic example of vertex calculation - a triangle
fig_num = 1002;
figure(fig_num);
clf;

fig_num = 2;
vertices = [0 0; 1 1; 0 1; 0 0]*5;
test_polytope.vertices = vertices;

% fill in other fields from the vertices field
test_polytope = fcn_VSkel_fillPolytopeFieldsFromVertices(test_polytope);

edge_cut = 0.1;
[shrunk_polytope, new_vertices, new_projection_vectors, cut_distance] = fcn_VSkel_polytopeShrinkFromEdges(...
    test_polytope,edge_cut,fig_num);

% Check variable types
assert(isstruct(shrunk_polytope));
assert(iscell(new_vertices));
assert(iscell(new_projection_vectors));
assert(isnumeric(cut_distance));

% Check variable sizes
num_nested = 2; % This is the number of nested figures within the vertex skeleton
assert(length(new_vertices)==num_nested);
assert(length(new_projection_vectors)==num_nested);
assert(isscalar(cut_distance(1,:)));
assert(length(cut_distance(:,1))==num_nested);

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Basic example of vertex calculation - a triangle with too big a cut
fig_num = 1003;
figure(fig_num);
clf;

vertices = [0 0; 1 1; 0 1; 0 0];
test_polytope.vertices = vertices;

% fill in other fields from the vertices field
test_polytope = fcn_VSkel_fillPolytopeFieldsFromVertices(test_polytope);

edge_cut = 2;
[shrunk_polytope, new_vertices, new_projection_vectors, cut_distance] = fcn_VSkel_polytopeShrinkFromEdges(...
    test_polytope,edge_cut,fig_num);


% Check variable types
assert(isstruct(shrunk_polytope));
assert(iscell(new_vertices));
assert(iscell(new_projection_vectors));
assert(isnumeric(cut_distance));

% Check variable sizes
num_nested = 2; % This is the number of nested figures within the vertex skeleton
assert(length(new_vertices)==num_nested);
assert(length(new_projection_vectors)==num_nested);
assert(isscalar(cut_distance(1,:)));
assert(length(cut_distance(:,1))==num_nested);

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));


%% Basic example of vertex calculation - a square, NO FIGURE
fig_num = 1004;
figure(fig_num);
close(fig_num);

vertices = [0 0; 1 0; 1 1; 0 1; 0 0]*5;
% vertices = [0 0; 2 0; 1 2; 0 1; 0 0]*5;
test_polytope.vertices = vertices;

% fill in other fields from the vertices field
test_polytope = fcn_VSkel_fillPolytopeFieldsFromVertices(test_polytope);


edge_cut = 0.1;
[shrunk_polytope, new_vertices, new_projection_vectors, cut_distance] = fcn_VSkel_polytopeShrinkFromEdges( test_polytope,edge_cut,[]);

% Check variable types
assert(isstruct(shrunk_polytope));
assert(iscell(new_vertices));
assert(iscell(new_projection_vectors));
assert(isnumeric(cut_distance));

% Check variable sizes
num_nested = 2; % This is the number of nested figures within the vertex skeleton
assert(length(new_vertices)==num_nested);
assert(length(new_projection_vectors)==num_nested);
assert(isscalar(cut_distance(1,:)));
assert(length(cut_distance(:,1))==num_nested);


% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));



%% Basic example of vertex calculation - a square, NO FIGURE, FAST MODE
fig_num = 1004;
figure(fig_num);
close(fig_num);

vertices = [0 0; 1 0; 1 1; 0 1; 0 0]*5;
% vertices = [0 0; 2 0; 1 2; 0 1; 0 0]*5;
test_polytope.vertices = vertices;

% fill in other fields from the vertices field
test_polytope = fcn_VSkel_fillPolytopeFieldsFromVertices(test_polytope);


edge_cut = 0.1;
[shrunk_polytope, new_vertices, new_projection_vectors, cut_distance] = fcn_VSkel_polytopeShrinkFromEdges( test_polytope,edge_cut,-1);

% Check variable types
assert(isstruct(shrunk_polytope));
assert(iscell(new_vertices));
assert(iscell(new_projection_vectors));
assert(isnumeric(cut_distance));

% Check variable sizes
num_nested = 2; % This is the number of nested figures within the vertex skeleton
assert(length(new_vertices)==num_nested);
assert(length(new_projection_vectors)==num_nested);
assert(isscalar(cut_distance(1,:)));
assert(length(cut_distance(:,1))==num_nested);


% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));


%% Compare speeds of pre-calculation versus post-calculation versus a fast variant
fig_num = 1006;
figure(fig_num);
clf;

figure(fig_num);
axis equal;
hold on;

% Set up polytopes
load('testData_fcn_VSkel_plotPolytopes.mat','polytopes','polytopes2');

% Pick a random polytope
Npolys = length(polytopes);
rand_poly = 1+floor(rand*Npolys);
shrinker = polytopes(rand_poly);

edge_cut = 0.1;
[shrunk_polytope, new_vertices, new_projection_vectors, cut_distance] = ...
    fcn_VSkel_polytopeShrinkFromEdges(...
        shrinker,edge_cut,fig_num);

edge_cut_step = 0.005;
iterations = (shrinker.max_radius/1.5+edge_cut_step)/edge_cut_step;

% Do calculation without pre-calculation
tic;
% for edge_cut = edge_cut_step:edge_cut_step:(shrinker.max_radius+edge_cut_step)
for edge_cut = edge_cut_step:edge_cut_step:(shrinker.max_radius/1.5+edge_cut_step)
    fcn_VSkel_polytopeShrinkFromEdges(...
        shrinker,edge_cut, -1);
end
slow_method = toc;

% Do calculation with pre-calculation, FAST_MODE on
tic;
% for edge_cut = edge_cut_step:edge_cut_step:(shrinker.max_radius+edge_cut_step)
for edge_cut = edge_cut_step:edge_cut_step:(shrinker.max_radius/1.5+edge_cut_step)
    fcn_VSkel_polytopeShrinkFromEdges(shrinker,edge_cut,new_vertices, new_projection_vectors, cut_distance, -1);
end
fast_method = toc;

% Do calculation with pre-calculation and ONLY vertex calculations,
% FAST_MODE on
tic;
% for edge_cut = edge_cut_step:edge_cut_step:(shrinker.max_radius+edge_cut_step)
for edge_cut = edge_cut_step:edge_cut_step:(shrinker.max_radius/1.5+edge_cut_step)
    fcn_VSkel_polytopeShrinkFromEdges_fast(shrinker,edge_cut,new_vertices, new_projection_vectors, cut_distance, -1);
end
fastest_method = toc;

% Plot results as bar chart
figure(373737);
clf;
X = categorical({'No Reuse','Reuse Skeleton','Fast reuse'});
X = reordercats(X,{'No Reuse','Reuse Skeleton','Fast reuse'});
Y = [slow_method fast_method fastest_method]*1000/iterations;
bar(X,Y)
ylabel('Execution time (Milliseconds)')