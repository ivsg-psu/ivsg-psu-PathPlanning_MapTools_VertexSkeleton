% script_test_fcn_VSkel_fcn_polytopeFindUnitDirectionVectors
% Tests function: fcn_VSkel_fcn_polytopeFindUnitDirectionVectors

% REVISION HISTORY:
% 2022_02_15
% -- first written by S. Brennan
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

%% Demonstration case 1: random 2D polytope
fig_num = 9001;
figure(fig_num);
clf;

vertices = [0 0; 2 0; 1 2; 0 1; 0 0]*5;
[unit_normal_vectors, vertex_projection_vectors] = fcn_VSkel_fcn_polytopeFindUnitDirectionVectors(vertices,fig_num);

% Check variable types
assert(length(unit_normal_vectors(:,1)) == length(vertices(:,1)));
assert(length(vertex_projection_vectors(:,1)) == length(vertices(:,1)));

% Check that all are unit vectors
unit_normal_vectors_length = sum(unit_normal_vectors.^2,2).^0.5;
vertex_projection_vectors_length = sum(unit_normal_vectors.^2,2).^0.5;

assert(all(abs(unit_normal_vectors_length - 1)<1E-10)==true);
assert(all(abs(vertex_projection_vectors_length - 1)<1E-10)==true);

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

%% Basic example of vertex calculation - non-normal wall shrinking
fig_num = 1001;
figure(fig_num);
clf;

% this polytope has a vertical wall
vertices = [0 0; 2 0; 1 2; 0 1; 0 0]*5;
[unit_normal_vectors, vertex_projection_vectors] = fcn_VSkel_fcn_polytopeFindUnitDirectionVectors(vertices,fig_num); %#ok<*ASGLU>

% Check variable types
assert(length(unit_normal_vectors(:,1)) == length(verticies(:,1)));
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


%% Example case 2: square - will have 2 nested solutions
fig_num = 1002;
figure(fig_num);
clf;

vertices = [0 0; 1 0; 1 1; 0 1; 0 0]*10;
[unit_normal_vectors, vertex_projection_vectors] = fcn_VSkel_fcn_polytopeFindUnitDirectionVectors(vertices,fig_num);


% assert that final vertices are within 5% error of having the same x
% position
final_vertices = new_vertices{end}(1,:);
expected_final_vertices = [5 5];
error_tolerance = 0.05;
vertical_error = sum(sum((final_vertices - expected_final_vertices).^2,1).^0.5,2);
% note this error is currently about 33% as of the commit this line was
% added
assert(vertical_error <= error_tolerance,['Expected and actual final points are: ',...
    '\n\t(%d,%d)  (expected) \n\t(%d,%d) (actual) \nyielding a vertical error of %d. Error tolerance was %d.\n'],...
    expected_final_vertices(1,1), expected_final_vertices(1,2),...
    final_vertices(1,1),final_vertices(1,2),vertical_error,error_tolerance);

% Check variable types
assert(length(unit_normal_vectors(:,1)) == length(verticies(:,1)));
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

%% Example case 3: wide rectangle
fig_num = 1003;
figure(fig_num);
clf;

vertices = [0 0; 1 0; 1 0.5; 0 0.5; 0 0]*10;
[unit_normal_vectors, vertex_projection_vectors] = fcn_VSkel_fcn_polytopeFindUnitDirectionVectors(vertices,fig_num);


% assert that final vertices are within 5% error of having the same x
% position
final_vertices = new_vertices{end}(1,:);
expected_final_vertices = [5 2.5];
error_tolerance = 0.05;
vertical_error = sum(sum((final_vertices - expected_final_vertices).^2,1).^0.5,2);
% note this error is currently about 33% as of the commit this line was
% added
assert(vertical_error <= error_tolerance,['Expected and actual final points are: ',...
    '\n\t(%d,%d)  (expected) \n\t(%d,%d) (actual) \nyielding a vertical error of %d. Error tolerance was %d.\n'],...
    expected_final_vertices(1,1), expected_final_vertices(1,2),...
    final_vertices(1,1),final_vertices(1,2),vertical_error,error_tolerance);


% Check variable types
assert(length(unit_normal_vectors(:,1)) == length(verticies(:,1)));
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

%% Example case 4: tall rectangle
fig_num = 1004;
figure(fig_num);
clf;

vertices = [0 0; 0.5 0; 0.5 1; 0 1; 0 0]*10;
[unit_normal_vectors, vertex_projection_vectors] = fcn_VSkel_fcn_polytopeFindUnitDirectionVectors(vertices,fig_num);

% assert that final vertices are within 5% error of having the same x
% position
final_vertices = new_vertices{end}(1,:);
expected_final_vertices = [2.5 5];
error_tolerance = 0.05;
vertical_error = sum(sum((final_vertices - expected_final_vertices).^2,1).^0.5,2);
% note this error is currently about 33% as of the commit this line was
% added
assert(vertical_error <= error_tolerance,['Expected and actual final points are: ',...
    '\n\t(%d,%d)  (expected) \n\t(%d,%d) (actual) \nyielding a vertical error of %d. Error tolerance was %d.\n'],...
    expected_final_vertices(1,1), expected_final_vertices(1,2),...
    final_vertices(1,1),final_vertices(1,2),vertical_error,error_tolerance);


% Check variable types
assert(length(unit_normal_vectors(:,1)) == length(verticies(:,1)));
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

%% Example case 5: goofy polytope
fig_num = 1005;
figure(fig_num);
clf;

vertices = [0 0; 10 0; 5 15; 4 17; 1 13; 0 5; 0 0];
[unit_normal_vectors, vertex_projection_vectors] = fcn_VSkel_fcn_polytopeFindUnitDirectionVectors(vertices,fig_num);

% Check variable types
assert(length(unit_normal_vectors(:,1)) == length(verticies(:,1)));
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

%% Random polytope calculation
fig_num = 1006;
figure(fig_num);
clf;

rng(1); % Set the random number generator

% Set up polytopes
load('testData_fcn_VSkel_plotPolytopes.mat','polytopes','polytopes2');

% Pick a random polytope
Npolys = length(polytopes);
rand_poly = 1+floor(rand*Npolys);
shrinker = polytopes(rand_poly);

% Do skeleton calculation
[unit_normal_vectors, vertex_projection_vectors] = fcn_VSkel_fcn_polytopeFindUnitDirectionVectors(shrinker.vertices,fig_num);

% Check variable types
assert(length(unit_normal_vectors(:,1)) == length(verticies(:,1)));
assert(iscell(new_projection_vectors));
assert(isnumeric(cut_distance));

% Check variable sizes
num_nested = length(new_vertices); % This is the number of nested figures within the vertex skeleton
assert(length(new_vertices)==num_nested);
assert(length(new_projection_vectors)==num_nested);
assert(isscalar(cut_distance(1,:)));
assert(length(cut_distance(:,1))==num_nested);

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

fig_num = 1106;
figure(fig_num); 
clf;
axis equal;
hold on;

% for edge_cut = edge_cut_step:edge_cut_step:(shrinker.max_radius+edge_cut_step)
edge_cut_step = 0.002;
for cut = edge_cut_step:edge_cut_step:(shrinker.max_radius+edge_cut_step)
    
    % Find the shape that is less than or equal to the cut
    shape_index = find(cut_distance<=cut,1,'last');
    
    % Grab vertices to start from, cut to start from
    starting_vertices = new_vertices{shape_index};
    starting_cut = cut_distance(shape_index);
    
    % Calculate projection distance
    projection_distance = cut - starting_cut;
    
    % Determine final vertices
    final_vertices = starting_vertices + new_projection_vectors{shape_index}*projection_distance;
    
    % Plot results
    plot(final_vertices(:,1),final_vertices(:,2),'.-','Linewidth',2,'Markersize',20);
end


%% Basic example of vertex calculation - non-normal wall shrinking, NO FIGURE
fig_num = 1007;
figure(fig_num);
close(fig_num);

% this polytope has a vertical wall
vertices = [0 0; 2 0; 1 2; 0 1; 0 0]*5;
[unit_normal_vectors, vertex_projection_vectors] = fcn_VSkel_fcn_polytopeFindUnitDirectionVectors(vertices,[]); %#ok<*ASGLU>

% Check variable types
assert(length(unit_normal_vectors(:,1)) == length(verticies(:,1)));
assert(iscell(new_projection_vectors));
assert(isnumeric(cut_distance));

% Check variable sizes
num_nested = 3; % This is the number of nested figures within the vertex skeleton
assert(length(new_vertices)==num_nested);
assert(length(new_projection_vectors)==num_nested);
assert(isscalar(cut_distance(1,:)));
assert(length(cut_distance(:,1))==num_nested);


% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));

%% Basic example of vertex calculation - non-normal wall shrinking, NO FIGURE, FAST MODE
fig_num = 1008;
figure(fig_num);
close(fig_num);

% this polytope has a vertical wall
vertices = [0 0; 2 0; 1 2; 0 1; 0 0]*5;
[unit_normal_vectors, vertex_projection_vectors] = fcn_VSkel_fcn_polytopeFindUnitDirectionVectors(vertices,-1); %#ok<*ASGLU>

% Check variable types
assert(length(unit_normal_vectors(:,1)) == length(verticies(:,1)));
assert(iscell(new_projection_vectors));
assert(isnumeric(cut_distance));

% Check variable sizes
num_nested = 3; % This is the number of nested figures within the vertex skeleton
assert(length(new_vertices)==num_nested);
assert(length(new_projection_vectors)==num_nested);
assert(isscalar(cut_distance(1,:)));
assert(length(cut_distance(:,1))==num_nested);


% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));