% script_test_fcn_VSkel_polytopeFindVertexSkeleton
% Tests function: fcn_VSkel_polytopeFindVertexSkeleton

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

%% Demonstration case 1: non-normal wall shrinking
fig_num = 9001;
figure(fig_num);
clf;

vertices = [0 0; 2 0; 1 2; 0 1; 0 0]*5;
[new_vertices, new_projection_vectors, cut_distance] = fcn_VSkel_polytopeFindVertexSkeleton(vertices,fig_num);

% Check variable types
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

% Make a plot of nested cuts?
if 1==0
    figure(4747);
    grid on
    grid minor
    hold on
    axis equal

    for cut = 0:0.5:cut_distance(end)

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
        plot(final_vertices(:,1),final_vertices(:,2),'r.-','Linewidth',2,'Markersize',20);
    end
end

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

%% Basic example of vertex calculation - three point test, horizontal line segment
fig_num = 1001;
figure(fig_num);
clf;

vertices = [0 0; 2 0; 0 0];
[new_vertices, new_projection_vectors, cut_distance] = fcn_VSkel_polytopeFindVertexSkeleton(vertices,fig_num); %#ok<*ASGLU>

% Check variable types
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

%% Basic example of vertex calculation - three point test, vertical line segment
fig_num = 1002;
figure(fig_num);
clf;

vertices = [0 0; 0 2; 0 0];
[new_vertices, new_projection_vectors, cut_distance] = fcn_VSkel_polytopeFindVertexSkeleton(vertices,fig_num); %#ok<*ASGLU>

% Check variable types
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

%% Basic example of vertex calculation - three point test, arbitrary line segment
fig_num = 1003;
figure(fig_num);
clf;

vertices = [-1 2; 3 4; -1 2];
[new_vertices, new_projection_vectors, cut_distance] = fcn_VSkel_polytopeFindVertexSkeleton(vertices,fig_num); %#ok<*ASGLU>

% Check variable types
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



%% Basic example of vertex calculation - non-normal wall shrinking

URHERE

fig_num = 1001;
figure(fig_num);
clf;

% this polytope has a vertical wall
vertices = [0 0; 2 0; 1 2; 0 1; 0 0]*5;
[new_vertices, new_projection_vectors, cut_distance] = fcn_VSkel_polytopeFindVertexSkeleton(vertices,fig_num); %#ok<*ASGLU>

% Check variable types
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


%% Example case 2: square - will have 2 nested solutions
fig_num = 1002;
figure(fig_num);
clf;

vertices = [0 0; 1 0; 1 1; 0 1; 0 0]*10;
[new_vertices, new_projection_vectors, cut_distance] = fcn_VSkel_polytopeFindVertexSkeleton(vertices,fig_num);


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

%% Example case 3: wide rectangle
fig_num = 1003;
figure(fig_num);
clf;

vertices = [0 0; 1 0; 1 0.5; 0 0.5; 0 0]*10;
[new_vertices, new_projection_vectors, cut_distance] = fcn_VSkel_polytopeFindVertexSkeleton(vertices,fig_num);


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

%% Example case 4: tall rectangle
fig_num = 1004;
figure(fig_num);
clf;

vertices = [0 0; 0.5 0; 0.5 1; 0 1; 0 0]*10;
[new_vertices, new_projection_vectors, cut_distance] = fcn_VSkel_polytopeFindVertexSkeleton(vertices,fig_num);

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

%% Example case 5: goofy polytope
fig_num = 1005;
figure(fig_num);
clf;

vertices = [0 0; 10 0; 5 15; 4 17; 1 13; 0 5; 0 0];
[new_vertices, new_projection_vectors, cut_distance] = fcn_VSkel_polytopeFindVertexSkeleton(vertices,fig_num);

% Check variable types
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
[new_vertices, new_projection_vectors, cut_distance] = fcn_VSkel_polytopeFindVertexSkeleton(shrinker.vertices,fig_num);

% Check variable types
assert(iscell(new_vertices));
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
[new_vertices, new_projection_vectors, cut_distance] = fcn_VSkel_polytopeFindVertexSkeleton(vertices,[]); %#ok<*ASGLU>

% Check variable types
assert(iscell(new_vertices));
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
[new_vertices, new_projection_vectors, cut_distance] = fcn_VSkel_polytopeFindVertexSkeleton(vertices,-1); %#ok<*ASGLU>

% Check variable types
assert(iscell(new_vertices));
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