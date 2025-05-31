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

%% Demonstration example 1: no merged edges
fig_num = 0001;
figure(fig_num);
clf;

clear vertices
vertices{1} = [0 0; 1 0; 1 1];
vertices{2} = [2 2; 3 4; 1 5];
polytopeStructure = fcn_VSkel_polytopeFillStructureFromVertices(vertices, (-1));
Nfaces = length(polytopeStructure.polyPatch.Faces(:,1));
randAngles = rand(Nfaces,1)*2*pi;
randomFaceVectors = [cos(randAngles) sin(randAngles)];
Nvertices = length(polytopeStructure.polyPatch.Vertices(:,1));
randAngles = rand(Nvertices,1)*2*pi;
randomVertexVectors = [cos(randAngles) sin(randAngles)];



clear plot_formatting
plot_formatting.vertices_plot.vertexLabels_flagOn = 1;
plot_formatting.vertices_plot.faceLabels_flagOn = 1;

h_fig =  fcn_VSkel_plotPolytopeDetails(...
       polytopeStructure,...
       (randomFaceVectors), ...  % unit_normal_vectors
       (randomVertexVectors), ...  % unit_vertex_projection_vectors
       (plot_formatting),... % plot_formatting
       (fig_num));  % fig_num

% Check variable types
assert(ishandle(h_fig));

% Check that it is the correct figure
assert(isequal(h_fig.Number,fig_num));


%% Demonstration case 2: plotting as in fcn_VSkel_polytopeFindUnitDirectionVectors
fig_num = 0002;
figure(fig_num);
clf;

vertices = [0 0; 2 0; 1 2; 0 1]*5;
polytopeStructure = fcn_VSkel_polytopeFillStructureFromVertices(vertices, (-1));
[unit_normal_vectors, vector_direction_of_unit_cut, ~] = ...
    fcn_VSkel_polytopeFindUnitDirectionVectors(polytopeStructure,-1);
clear plot_formatting
plot_formatting.vertices_plot.vertexLabels_flagOn = 1;
plot_formatting.vertices_plot.faceLabels_flagOn = 1;

h_fig =  fcn_VSkel_plotPolytopeDetails(...
       polytopeStructure,...
       (unit_normal_vectors), ...  % unit_normal_vectors
       (vector_direction_of_unit_cut), ... % vector_direction_of_unit_cut
       (plot_formatting),... % plot_formatting
       (fig_num));  % fig_num


% Check variable types
assert(ishandle(h_fig));

% Check that it is the correct figure
assert(isequal(h_fig.Number,fig_num));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Demonstration case 3: merged edges in 2D
fig_num = 0003;
figure(fig_num);
clf;

clear vertices
vertices{1} = [0 0; 1 0; 1 1];
vertices{2} = [0 0; 1 1; 0 1];
polytopeStructure = fcn_VSkel_polytopeFillStructureFromVertices(vertices, (-1));
Nfaces = length(polytopeStructure.polyPatch.Faces(:,1));
randAngles = rand(Nfaces,1)*2*pi;
randomFaceVectors = [cos(randAngles) sin(randAngles)];
Nvertices = length(polytopeStructure.polyPatch.Vertices(:,1));
randAngles = rand(Nvertices,1)*2*pi;
randomVertexVectors = [cos(randAngles) sin(randAngles)];



clear plot_formatting
plot_formatting.vertices_plot.vertexLabels_flagOn = 1;
plot_formatting.vertices_plot.faceLabels_flagOn = 1;

h_fig =  fcn_VSkel_plotPolytopeDetails(...
       polytopeStructure,...
       (randomFaceVectors), ...  % unit_normal_vectors
       (randomVertexVectors), ...  % unit_vertex_projection_vectors
       (plot_formatting),... % plot_formatting
       (fig_num));  % fig_num

% Check variable types
assert(ishandle(h_fig));

% Check that it is the correct figure
assert(isequal(h_fig.Number,fig_num));

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

vertices = [0 0; 2 0; 1 2; 0 1]*5;
polytopeStructure = fcn_VSkel_polytopeFillStructureFromVertices(vertices, (-1));
[unit_normal_vectors, vector_direction_of_unit_cut, ~] = ...
    fcn_VSkel_polytopeFindUnitDirectionVectors(polytopeStructure,-1);
clear plot_formatting
plot_formatting.vertices_plot.vertexLabels_flagOn = 1;
plot_formatting.vertices_plot.faceLabels_flagOn = 1;


plot_formatting.vertices_plot.faceLabels_Color = [0 0 0]; % Change text to black


h_fig =  fcn_VSkel_plotPolytopeDetails(...
       polytopeStructure,...
       (unit_normal_vectors), ...  % unit_normal_vectors
       (vector_direction_of_unit_cut), ... % vector_direction_of_unit_cut
       (plot_formatting),... % plot_formatting
       (fig_num));  % fig_num


% Check variable types
assert(ishandle(h_fig));

% Check that it is the correct figure
assert(isequal(h_fig.Number,fig_num));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));


%% Demonstration case 2: one point
fig_num = 1002;
figure(fig_num);
clf;

vertices = [0 2; 0 2]; % A single point
polytopeStructure = fcn_VSkel_polytopeFillStructureFromVertices(vertices, (-1));
[unit_normal_vectors, vector_direction_of_unit_cut, ~] = ...
    fcn_VSkel_polytopeFindUnitDirectionVectors(polytopeStructure,-1);
clear plot_formatting
plot_formatting.vertices_plot.vertexLabels_flagOn = 1;
plot_formatting.vertices_plot.faceLabels_flagOn = 1;


plot_formatting.vertices_plot.faceLabels_Color = [0 0 0]; % Change text to black


h_fig =  fcn_VSkel_plotPolytopeDetails(...
       polytopeStructure,...
       (unit_normal_vectors), ...  % unit_normal_vectors
       (vector_direction_of_unit_cut), ... % vector_direction_of_unit_cut
       (plot_formatting),... % plot_formatting
       (fig_num));  % fig_num


% Check variable types
assert(ishandle(h_fig));

% Check that it is the correct figure
assert(isequal(h_fig.Number,fig_num));

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
