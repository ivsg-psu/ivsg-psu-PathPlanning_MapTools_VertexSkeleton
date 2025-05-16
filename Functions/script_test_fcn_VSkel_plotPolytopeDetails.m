% script_test_fcn_VSkel_plotPolytopeDetails
% Tests function: fcn_VSkel_plotPolytopeDetails

% REVISION HISTORY:
% 2025_05_04 by Sean Brennan
% -- first written by S. Brennan
% 2025_05_15 by Sean Brennan
% -- added case where vertices can be only one point

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

%% Demonstration case 1: plotting as in fcn_VSkel_polytopeFindUnitDirectionVectors
fig_num = 0001;
figure(fig_num);
clf;

vertices = [0 0; 2 0; 1 2; 0 1; 0 0]*5;
[unit_normal_vectors, unit_vertex_projection_vectors, vector_direction_of_unit_cut, flag_vertexIsNonConvex] = ...
    fcn_VSkel_polytopeFindUnitDirectionVectors(vertices,-1);

h_fig =  fcn_VSkel_plotPolytopeDetails(...
       vertices,...
       (unit_normal_vectors), ...  % unit_normal_vectors
       (unit_vertex_projection_vectors), ...  % unit_vertex_projection_vectors
       (vector_direction_of_unit_cut), ... % vector_direction_of_unit_cut
       (flag_vertexIsNonConvex),...  % flag_vertexIsNonConvex
       (1),...  % flag_plotEdgeGhostlines
       (1),...  % flag_plotVertexProjectionGhostlines
       ([]),... % plot_formatting
       (fig_num));  % fig_num


% Check variable types
assert(length(unit_normal_vectors(:,1)) == length(vertices(:,1)));
assert(length(unit_vertex_projection_vectors(:,1)) == length(vertices(:,1)));
assert(length(vector_direction_of_unit_cut(:,1)) == length(vertices(:,1)));
assert(length(flag_vertexIsNonConvex(:,1)) == length(vertices(:,1)));

% Check that all unit vectors are unit length
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

%% Setting plot styles
fig_num = 1001;
figure(fig_num);
clf;

vertices = [0 0; 2 0; 1 2; 0 1; 0 0]*5;
[unit_normal_vectors, unit_vertex_projection_vectors, vector_direction_of_unit_cut, flag_vertexIsNonConvex] = ...
    fcn_VSkel_polytopeFindUnitDirectionVectors(vertices,-1);

plot_formatting.vertices_plot.edgeLabelsColor = [0 0 0]; % Change text to black

h_fig =  fcn_VSkel_plotPolytopeDetails(...
       vertices,...
       (unit_normal_vectors), ...  % unit_normal_vectors
       (unit_vertex_projection_vectors), ...  % unit_vertex_projection_vectors
       (vector_direction_of_unit_cut), ... % vector_direction_of_unit_cut
       (flag_vertexIsNonConvex),...  % flag_vertexIsNonConvex
       (1),...  % flag_plotEdgeGhostlines
       (1),...  % flag_plotVertexProjectionGhostlines
       (plot_formatting),...  % plot_formatting
       (fig_num));  % fig_num

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Demonstration case 2: one point
fig_num = 0002;
figure(fig_num);
clf;

vertices = [0 2; 0 2]; % A single point

[unit_normal_vectors, unit_vertex_projection_vectors, vector_direction_of_unit_cut, flag_vertexIsNonConvex] = ...
    fcn_VSkel_polytopeFindUnitDirectionVectors(vertices,-1);

h_fig =  fcn_VSkel_plotPolytopeDetails(...
       vertices,...
       (unit_normal_vectors), ...  % unit_normal_vectors
       (unit_vertex_projection_vectors), ...  % unit_vertex_projection_vectors
       (vector_direction_of_unit_cut), ... % vector_direction_of_unit_cut
       (flag_vertexIsNonConvex),...  % flag_vertexIsNonConvex
       (1),...  % flag_plotEdgeGhostlines
       (1),...  % flag_plotVertexProjectionGhostlines
       ([]),... % plot_formatting
       (fig_num));  % fig_num


% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));


%% Fast Mode Tests
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  ______        _     __  __           _        _______        _
% |  ____|      | |   |  \/  |         | |      |__   __|      | |
% | |__ __ _ ___| |_  | \  / | ___   __| | ___     | | ___  ___| |_ ___
% |  __/ _` / __| __| | |\/| |/ _ \ / _` |/ _ \    | |/ _ \/ __| __/ __|
% | | | (_| \__ \ |_  | |  | | (_) | (_| |  __/    | |  __/\__ \ |_\__ \
% |_|  \__,_|___/\__| |_|  |_|\___/ \__,_|\___|    |_|\___||___/\__|___/
%
%
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Fast%20Mode%20Tests
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% None - cannot do fast mode on plotting routines
