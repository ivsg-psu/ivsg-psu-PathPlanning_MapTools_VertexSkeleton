% script_test_fcn_VSkel_polytopeFindMinimumEnclosedSphere
% Tests function: fcn_VSkel_polytopeFindMinimumEnclosedSphere

% REVISION HISTORY:
% 2025_05_03 by Sean Brennan
% -- first written by S. Brennan

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
fig_num = 0001;
figure(fig_num);
clf;

vertices = [0 0; 2 0; 1 2; 0 1; 0 0]*5;
[unit_normal_vectors, unit_vertex_projection_vectors, vector_direction_of_unit_cut, flag_vertexIsNonConvex]  = fcn_VSkel_polytopeFindUnitDirectionVectors(vertices,-1);
max_edge_cuts = fcn_VSkel_polytopeFindMaxEdgeCut(vertices, unit_normal_vectors, unit_vertex_projection_vectors, (-1));
[sphereRadii, definingBoundaries] = fcn_VSkel_polytopeFindEnclosedSpheres(vertices, unit_normal_vectors, unit_vertex_projection_vectors, vector_direction_of_unit_cut, flag_vertexIsNonConvex, max_edge_cuts, (-1));

[min_cut, boundaryEngagedAtMinCut, indices_repeated, intersection_points] = ...
    fcn_VSkel_polytopeFindMinimumEnclosedSphere(vertices, ...
    sphereRadii, definingBoundaries, unit_normal_vectors, unit_vertex_projection_vectors, ...
    vector_direction_of_unit_cut, (fig_num));

% Check variable types
assert(isscalar(min_cut));
assert(isnumeric(min_cut));
assert(isnumeric(boundaryEngagedAtMinCut));
assert(isnumeric(indices_repeated));
assert(isnumeric(intersection_points));

% Check lengths
assert(isequal(size(min_cut),[1 1]));
assert(length(boundaryEngagedAtMinCut) == length(indices_repeated(:,1)));
assert(length(indices_repeated(1,:))==1);
assert(length(indices_repeated(:,1))>=1);
assert(length(intersection_points(1,:))==2);
assert(length(intersection_points(:,1))==length(vertices(:,1)));


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

%% Demonstration case 1: random 2D polytope
fig_num = 1001;
figure(fig_num);
clf;

vertices = [0 0; 2 0; 1 2; 0 1; 0 0]*5;
[unit_normal_vectors, unit_vertex_projection_vectors, vector_direction_of_unit_cut, flag_vertexIsNonConvex]  = fcn_VSkel_polytopeFindUnitDirectionVectors(vertices,-1);
max_edge_cuts = fcn_VSkel_polytopeFindMaxEdgeCut(vertices, unit_normal_vectors, unit_vertex_projection_vectors, (-1));
[sphereRadii, definingBoundaries] = fcn_VSkel_polytopeFindEnclosedSpheres(vertices, unit_normal_vectors, unit_vertex_projection_vectors, vector_direction_of_unit_cut, flag_vertexIsNonConvex, max_edge_cuts, (-1));


[min_cut, boundaryEngagedAtMinCut, indices_repeated, intersection_points] = ...
    fcn_VSkel_polytopeFindMinimumEnclosedSphere(vertices, ...
    sphereRadii, definingBoundaries, unit_normal_vectors, unit_vertex_projection_vectors, ...
    vector_direction_of_unit_cut, (fig_num));

% Check variable types
assert(isscalar(min_cut));
assert(isnumeric(min_cut));
assert(isnumeric(boundaryEngagedAtMinCut));
assert(isnumeric(indices_repeated));
assert(isnumeric(intersection_points));

% Check lengths
assert(isequal(size(min_cut),[1 1]));
assert(length(boundaryEngagedAtMinCut) == length(indices_repeated(:,1)));
assert(length(indices_repeated(1,:))==1);
assert(length(indices_repeated(:,1))>=1);
assert(length(intersection_points(1,:))==2);
assert(length(intersection_points(:,1))==length(vertices(:,1)));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Example case 2: square
fig_num = 1002;
figure(fig_num);
clf;

vertices = [0 0; 1 0; 1 1; 0 1; 0 0]*10;
[unit_normal_vectors, unit_vertex_projection_vectors, vector_direction_of_unit_cut, flag_vertexIsNonConvex]  = fcn_VSkel_polytopeFindUnitDirectionVectors(vertices,-1);
max_edge_cuts = fcn_VSkel_polytopeFindMaxEdgeCut(vertices, unit_normal_vectors, unit_vertex_projection_vectors, (-1));
[sphereRadii, definingBoundaries] = fcn_VSkel_polytopeFindEnclosedSpheres(vertices, unit_normal_vectors, unit_vertex_projection_vectors, vector_direction_of_unit_cut, flag_vertexIsNonConvex, max_edge_cuts, (-1));

[min_cut, boundaryEngagedAtMinCut, indices_repeated, intersection_points] = ...
    fcn_VSkel_polytopeFindMinimumEnclosedSphere(vertices, ...
    sphereRadii, definingBoundaries, unit_normal_vectors, unit_vertex_projection_vectors, ...
    vector_direction_of_unit_cut, (fig_num));

% Check variable types
assert(isscalar(min_cut));
assert(isnumeric(min_cut));
assert(isnumeric(boundaryEngagedAtMinCut));
assert(isnumeric(indices_repeated));
assert(isnumeric(intersection_points));

% Check lengths
assert(isequal(size(min_cut),[1 1]));
assert(length(boundaryEngagedAtMinCut) == length(indices_repeated(:,1)));
assert(length(indices_repeated(1,:))==1);
assert(length(indices_repeated(:,1))>=1);
assert(length(intersection_points(1,:))==2);
assert(length(intersection_points(:,1))==length(vertices(:,1)));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));


%% Example case 3: wide rectangle
fig_num = 1003;
figure(fig_num);
clf;

vertices = [0 0; 1 0; 1 0.5; 0 0.5; 0 0]*10;
[unit_normal_vectors, unit_vertex_projection_vectors, vector_direction_of_unit_cut, flag_vertexIsNonConvex]  = fcn_VSkel_polytopeFindUnitDirectionVectors(vertices,-1);
max_edge_cuts = fcn_VSkel_polytopeFindMaxEdgeCut(vertices, unit_normal_vectors, unit_vertex_projection_vectors, (-1));
[sphereRadii, definingBoundaries] = fcn_VSkel_polytopeFindEnclosedSpheres(vertices, unit_normal_vectors, unit_vertex_projection_vectors, vector_direction_of_unit_cut, flag_vertexIsNonConvex, max_edge_cuts, (-1));

[min_cut, boundaryEngagedAtMinCut, indices_repeated, intersection_points] = ...
    fcn_VSkel_polytopeFindMinimumEnclosedSphere(vertices, ...
    sphereRadii, definingBoundaries, unit_normal_vectors, unit_vertex_projection_vectors, ...
    vector_direction_of_unit_cut, (fig_num));

% Check variable types
assert(isscalar(min_cut));
assert(isnumeric(min_cut));
assert(isnumeric(boundaryEngagedAtMinCut));
assert(isnumeric(indices_repeated));
assert(isnumeric(intersection_points));

% Check lengths
assert(isequal(size(min_cut),[1 1]));
assert(length(boundaryEngagedAtMinCut) == length(indices_repeated(:,1)));
assert(length(indices_repeated(1,:))==1);
assert(length(indices_repeated(:,1))>=1);
assert(length(intersection_points(1,:))==2);
assert(length(intersection_points(:,1))==length(vertices(:,1)));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Example case 4: tall rectangle
fig_num = 1004;
figure(fig_num);
clf;

vertices = [0 0; 0.5 0; 0.5 1; 0 1; 0 0]*10;
[unit_normal_vectors, unit_vertex_projection_vectors, vector_direction_of_unit_cut, flag_vertexIsNonConvex]  = fcn_VSkel_polytopeFindUnitDirectionVectors(vertices,-1);
max_edge_cuts = fcn_VSkel_polytopeFindMaxEdgeCut(vertices, unit_normal_vectors, unit_vertex_projection_vectors, (-1));
[sphereRadii, definingBoundaries] = fcn_VSkel_polytopeFindEnclosedSpheres(vertices, unit_normal_vectors, unit_vertex_projection_vectors, vector_direction_of_unit_cut, flag_vertexIsNonConvex, max_edge_cuts, (-1));

[min_cut, boundaryEngagedAtMinCut, indices_repeated, intersection_points] = ...
    fcn_VSkel_polytopeFindMinimumEnclosedSphere(vertices, ...
    sphereRadii, definingBoundaries, unit_normal_vectors, unit_vertex_projection_vectors, ...
    vector_direction_of_unit_cut, (fig_num));

% Check variable types
assert(isscalar(min_cut));
assert(isnumeric(min_cut));
assert(isnumeric(boundaryEngagedAtMinCut));
assert(isnumeric(indices_repeated));
assert(isnumeric(intersection_points));

% Check lengths
assert(isequal(size(min_cut),[1 1]));
assert(length(boundaryEngagedAtMinCut) == length(indices_repeated(:,1)));
assert(length(indices_repeated(1,:))==1);
assert(length(indices_repeated(:,1))>=1);
assert(length(intersection_points(1,:))==2);
assert(length(intersection_points(:,1))==length(vertices(:,1)));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Example case 5: goofy polytope
fig_num = 1005;
figure(fig_num);
clf;

vertices = [0 0; 10 0; 5 15; 4 17; 1 13; 0 5; 0 0];
[unit_normal_vectors, unit_vertex_projection_vectors, vector_direction_of_unit_cut, flag_vertexIsNonConvex]  = fcn_VSkel_polytopeFindUnitDirectionVectors(vertices,-1);
max_edge_cuts = fcn_VSkel_polytopeFindMaxEdgeCut(vertices, unit_normal_vectors, unit_vertex_projection_vectors, (-1));
[sphereRadii, definingBoundaries] = fcn_VSkel_polytopeFindEnclosedSpheres(vertices, unit_normal_vectors, unit_vertex_projection_vectors, vector_direction_of_unit_cut, flag_vertexIsNonConvex, max_edge_cuts, (-1));

[min_cut, boundaryEngagedAtMinCut, indices_repeated, intersection_points] = ...
    fcn_VSkel_polytopeFindMinimumEnclosedSphere(vertices, ...
    sphereRadii, definingBoundaries, unit_normal_vectors, unit_vertex_projection_vectors, ...
    vector_direction_of_unit_cut, (fig_num));

% Check variable types
assert(isscalar(min_cut));
assert(isnumeric(min_cut));
assert(isnumeric(boundaryEngagedAtMinCut));
assert(isnumeric(indices_repeated));
assert(isnumeric(intersection_points));

% Check lengths
assert(isequal(size(min_cut),[1 1]));
assert(length(boundaryEngagedAtMinCut) == length(indices_repeated(:,1)));
assert(length(indices_repeated(1,:))==1);
assert(length(indices_repeated(:,1))>=1);
assert(length(intersection_points(1,:))==2);
assert(length(intersection_points(:,1))==length(vertices(:,1)));


% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Basic example of vertex calculation - non-convex 2D polytope
fig_num = 1006;
figure(fig_num);
clf;

% this polytope has a vertical wall
vertices = [0 0; 3/5 0; 1 1; 7/5 0; 2 0; 1 2; 0 1; 0 0]*5;
[unit_normal_vectors, unit_vertex_projection_vectors, vector_direction_of_unit_cut, flag_vertexIsNonConvex]  = fcn_VSkel_polytopeFindUnitDirectionVectors(vertices,-1);
max_edge_cuts = fcn_VSkel_polytopeFindMaxEdgeCut(vertices, unit_normal_vectors, unit_vertex_projection_vectors, (-1));
[sphereRadii, definingBoundaries] = fcn_VSkel_polytopeFindEnclosedSpheres(vertices, unit_normal_vectors, unit_vertex_projection_vectors, vector_direction_of_unit_cut, flag_vertexIsNonConvex, max_edge_cuts, (-1));

[min_cut, boundaryEngagedAtMinCut, indices_repeated, intersection_points] = ...
    fcn_VSkel_polytopeFindMinimumEnclosedSphere(vertices, ...
    sphereRadii, definingBoundaries, unit_normal_vectors, unit_vertex_projection_vectors, ...
    vector_direction_of_unit_cut, (fig_num));

% Check variable types
assert(isscalar(min_cut));
assert(isnumeric(min_cut));
assert(isnumeric(boundaryEngagedAtMinCut));
assert(isnumeric(indices_repeated));
assert(isnumeric(intersection_points));

% Check lengths
assert(isequal(size(min_cut),[1 1]));
assert(length(boundaryEngagedAtMinCut) == length(indices_repeated(:,1)));
assert(length(indices_repeated(1,:))==1);
assert(length(indices_repeated(:,1))>=1);
assert(length(intersection_points(1,:))==2);
assert(length(intersection_points(:,1))==length(vertices(:,1)));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Basic example of vertex calculation - hard non-convex 2D polytope
fig_num = 1007;
figure(fig_num);
clf;

% this polytope has a vertical wall
vertices = [0 0; 5 0; 6 4; 7 0; 10 0; 10 10; 4 10; 4 5; 0 5; 0 0];
[unit_normal_vectors, unit_vertex_projection_vectors, vector_direction_of_unit_cut, flag_vertexIsNonConvex]  = fcn_VSkel_polytopeFindUnitDirectionVectors(vertices,-1);
max_edge_cuts = fcn_VSkel_polytopeFindMaxEdgeCut(vertices, unit_normal_vectors, unit_vertex_projection_vectors, (-1));
[sphereRadii, definingBoundaries] = fcn_VSkel_polytopeFindEnclosedSpheres(vertices, unit_normal_vectors, unit_vertex_projection_vectors, vector_direction_of_unit_cut, flag_vertexIsNonConvex, max_edge_cuts, (-1));

[min_cut, boundaryEngagedAtMinCut, indices_repeated, intersection_points] = ...
    fcn_VSkel_polytopeFindMinimumEnclosedSphere(vertices, ...
    sphereRadii, definingBoundaries, unit_normal_vectors, unit_vertex_projection_vectors, ...
    vector_direction_of_unit_cut, (fig_num));

% Check variable types
assert(isscalar(min_cut));
assert(isnumeric(min_cut));
assert(isnumeric(boundaryEngagedAtMinCut));
assert(isnumeric(indices_repeated));
assert(isnumeric(intersection_points));

% Check lengths
assert(isequal(size(min_cut),[1 1]));
assert(length(boundaryEngagedAtMinCut) == length(indices_repeated(:,1)));
assert(length(indices_repeated(1,:))==1);
assert(length(indices_repeated(:,1))>=1);
assert(length(intersection_points(1,:))==2);
assert(length(intersection_points(:,1))==length(vertices(:,1)));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));


%% Basic example of vertex calculation - symmetric non-convex vertices
fig_num = 1008;
figure(fig_num);
clf;

% this polytope has two nonconvex vertices facing each other
vertices = [5 4; 6 0; 10 0; 10 10; 6 10; 5 6; 4 10; 0 10; 0 0; 4 0; 5 4];
[unit_normal_vectors, unit_vertex_projection_vectors, vector_direction_of_unit_cut, flag_vertexIsNonConvex]  = fcn_VSkel_polytopeFindUnitDirectionVectors(vertices,-1);
max_edge_cuts = fcn_VSkel_polytopeFindMaxEdgeCut(vertices, unit_normal_vectors, unit_vertex_projection_vectors, (-1));
[sphereRadii, definingBoundaries] = fcn_VSkel_polytopeFindEnclosedSpheres(vertices, unit_normal_vectors, unit_vertex_projection_vectors, vector_direction_of_unit_cut, flag_vertexIsNonConvex, max_edge_cuts, (-1));

[min_cut, boundaryEngagedAtMinCut, indices_repeated, intersection_points] = ...
    fcn_VSkel_polytopeFindMinimumEnclosedSphere(vertices, ...
    sphereRadii, definingBoundaries, unit_normal_vectors, unit_vertex_projection_vectors, ...
    vector_direction_of_unit_cut, (fig_num));

% Check variable types
assert(isscalar(min_cut));
assert(isnumeric(min_cut));
assert(isnumeric(boundaryEngagedAtMinCut));
assert(isnumeric(indices_repeated));
assert(isnumeric(intersection_points));

% Check lengths
assert(isequal(size(min_cut),[1 1]));
assert(length(boundaryEngagedAtMinCut) == length(indices_repeated(:,1)));
assert(length(indices_repeated(1,:))==1);
assert(length(indices_repeated(:,1))>=1);
assert(length(intersection_points(1,:))==2);
assert(length(intersection_points(:,1))==length(vertices(:,1)));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Basic example of vertex calculation - non-symmetric non-convex vertices
fig_num = 1009;
figure(fig_num);
clf;

% this polytope has two nonconvex vertices facing each other
vertices = [5 4; 6 0; 10 0; 10 7; 8 8; 10 9; 10 15; 0 15; 0 0; 4 0; 5 4];
[unit_normal_vectors, unit_vertex_projection_vectors, vector_direction_of_unit_cut, flag_vertexIsNonConvex]  = fcn_VSkel_polytopeFindUnitDirectionVectors(vertices,-1);
max_edge_cuts = fcn_VSkel_polytopeFindMaxEdgeCut(vertices, unit_normal_vectors, unit_vertex_projection_vectors, (-1));
[sphereRadii, definingBoundaries] = fcn_VSkel_polytopeFindEnclosedSpheres(vertices, unit_normal_vectors, unit_vertex_projection_vectors, vector_direction_of_unit_cut, flag_vertexIsNonConvex, max_edge_cuts, (-1));

[min_cut, boundaryEngagedAtMinCut, indices_repeated, intersection_points] = ...
    fcn_VSkel_polytopeFindMinimumEnclosedSphere(vertices, ...
    sphereRadii, definingBoundaries, unit_normal_vectors, unit_vertex_projection_vectors, ...
    vector_direction_of_unit_cut, (fig_num));

% Check variable types
assert(isscalar(min_cut));
assert(isnumeric(min_cut));
assert(isnumeric(boundaryEngagedAtMinCut));
assert(isnumeric(indices_repeated));
assert(isnumeric(intersection_points));

% Check lengths
assert(isequal(size(min_cut),[1 1]));
assert(length(boundaryEngagedAtMinCut) == length(indices_repeated(:,1)));
assert(length(indices_repeated(1,:))==1);
assert(length(indices_repeated(:,1))>=1);
assert(length(intersection_points(1,:))==2);
assert(length(intersection_points(:,1))==length(vertices(:,1)));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Basic example of vertex calculation - non-symmetric non-convex vertices
fig_num = 1010;
figure(fig_num);
clf;

% this polytope has two nonconvex vertices facing each other
vertices = [5 4; 6 0; 10 0; 10 7; 8 8; 10 9; 15 15; 0 15; 0 0; 4 0; 5 4];
[unit_normal_vectors, unit_vertex_projection_vectors, vector_direction_of_unit_cut, flag_vertexIsNonConvex]  = fcn_VSkel_polytopeFindUnitDirectionVectors(vertices,-1);
max_edge_cuts = fcn_VSkel_polytopeFindMaxEdgeCut(vertices, unit_normal_vectors, unit_vertex_projection_vectors, (-1));
[sphereRadii, definingBoundaries] = fcn_VSkel_polytopeFindEnclosedSpheres(vertices, unit_normal_vectors, unit_vertex_projection_vectors, vector_direction_of_unit_cut, flag_vertexIsNonConvex, max_edge_cuts, (-1));

[min_cut, boundaryEngagedAtMinCut, indices_repeated, intersection_points] = ...
    fcn_VSkel_polytopeFindMinimumEnclosedSphere(vertices, ...
    sphereRadii, definingBoundaries, unit_normal_vectors, unit_vertex_projection_vectors, ...
    vector_direction_of_unit_cut, (fig_num));

% Check variable types
assert(isscalar(min_cut));
assert(isnumeric(min_cut));
assert(isnumeric(boundaryEngagedAtMinCut));
assert(isnumeric(indices_repeated));
assert(isnumeric(intersection_points));

% Check lengths
assert(isequal(size(min_cut),[1 1]));
assert(length(boundaryEngagedAtMinCut) == length(indices_repeated(:,1)));
assert(length(indices_repeated(1,:))==1);
assert(length(indices_repeated(:,1))>=1);
assert(length(intersection_points(1,:))==2);
assert(length(intersection_points(:,1))==length(vertices(:,1)));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));


%% Basic example of vertex calculation - simple non-convex 2D polytope
fig_num = 1011;
figure(fig_num);
clf;

% this polytope has a vertical wall
vertices = [0 0; 3 0; 5 5; 7 0; 10 0; 10 6; 0 6; 0 0];
[unit_normal_vectors, unit_vertex_projection_vectors, vector_direction_of_unit_cut, flag_vertexIsNonConvex]  = fcn_VSkel_polytopeFindUnitDirectionVectors(vertices,-1);
max_edge_cuts = fcn_VSkel_polytopeFindMaxEdgeCut(vertices, unit_normal_vectors, unit_vertex_projection_vectors, (-1));
[sphereRadii, definingBoundaries] = fcn_VSkel_polytopeFindEnclosedSpheres(vertices, unit_normal_vectors, unit_vertex_projection_vectors, vector_direction_of_unit_cut, flag_vertexIsNonConvex, max_edge_cuts, (-1));

[min_cut, boundaryEngagedAtMinCut, indices_repeated, intersection_points] = ...
    fcn_VSkel_polytopeFindMinimumEnclosedSphere(vertices, ...
    sphereRadii, definingBoundaries, unit_normal_vectors, unit_vertex_projection_vectors, ...
    vector_direction_of_unit_cut, (fig_num));

% Check variable types
assert(isscalar(min_cut));
assert(isnumeric(min_cut));
assert(isnumeric(boundaryEngagedAtMinCut));
assert(isnumeric(indices_repeated));
assert(isnumeric(intersection_points));

% Check lengths
assert(isequal(size(min_cut),[1 1]));
assert(length(boundaryEngagedAtMinCut) == length(indices_repeated(:,1)));
assert(length(indices_repeated(1,:))==1);
assert(length(indices_repeated(:,1))>=1);
assert(length(intersection_points(1,:))==2);
assert(length(intersection_points(:,1))==length(vertices(:,1)));

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

%% Basic example of vertex calculation - non-normal wall shrinking, NO FIGURE
fig_num = 9901;
figure(fig_num);
close(fig_num);

vertices = [0 0; 2 0; 1 2; 0 1; 0 0]*5;
[unit_normal_vectors, unit_vertex_projection_vectors, vector_direction_of_unit_cut, flag_vertexIsNonConvex]  = fcn_VSkel_polytopeFindUnitDirectionVectors(vertices,-1);
max_edge_cuts = fcn_VSkel_polytopeFindMaxEdgeCut(vertices, unit_normal_vectors, unit_vertex_projection_vectors, (-1));
[sphereRadii, definingBoundaries] = fcn_VSkel_polytopeFindEnclosedSpheres(vertices, unit_normal_vectors, unit_vertex_projection_vectors, vector_direction_of_unit_cut, flag_vertexIsNonConvex, max_edge_cuts, (-1));

[min_cut, boundaryEngagedAtMinCut, indices_repeated, intersection_points] = ...
    fcn_VSkel_polytopeFindMinimumEnclosedSphere(vertices, ...
    sphereRadii, definingBoundaries, unit_normal_vectors, unit_vertex_projection_vectors, ...
    vector_direction_of_unit_cut, ([]));

% Check variable types
assert(isscalar(min_cut));
assert(isnumeric(min_cut));
assert(isnumeric(boundaryEngagedAtMinCut));
assert(isnumeric(indices_repeated));
assert(isnumeric(intersection_points));

% Check lengths
assert(isequal(size(min_cut),[1 1]));
assert(length(boundaryEngagedAtMinCut) == length(indices_repeated(:,1)));
assert(length(indices_repeated(1,:))==1);
assert(length(indices_repeated(:,1))>=1);
assert(length(intersection_points(1,:))==2);
assert(length(intersection_points(:,1))==length(vertices(:,1)));

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));

%% Basic example of vertex calculation - non-normal wall shrinking, NO FIGURE, FAST MODE
fig_num = 9902;
figure(fig_num);
close(fig_num);

vertices = [0 0; 2 0; 1 2; 0 1; 0 0]*5;
[unit_normal_vectors, unit_vertex_projection_vectors, vector_direction_of_unit_cut, flag_vertexIsNonConvex]  = fcn_VSkel_polytopeFindUnitDirectionVectors(vertices,-1);
max_edge_cuts = fcn_VSkel_polytopeFindMaxEdgeCut(vertices, unit_normal_vectors, unit_vertex_projection_vectors, (-1));
[sphereRadii, definingBoundaries] = fcn_VSkel_polytopeFindEnclosedSpheres(vertices, unit_normal_vectors, unit_vertex_projection_vectors, vector_direction_of_unit_cut, flag_vertexIsNonConvex, max_edge_cuts, (-1));

[min_cut, boundaryEngagedAtMinCut, indices_repeated, intersection_points] = ...
    fcn_VSkel_polytopeFindMinimumEnclosedSphere(vertices, ...
    sphereRadii, definingBoundaries, unit_normal_vectors, unit_vertex_projection_vectors, ...
    vector_direction_of_unit_cut, (-1));

% Check variable types
assert(isscalar(min_cut));
assert(isnumeric(min_cut));
assert(isnumeric(boundaryEngagedAtMinCut));
assert(isnumeric(indices_repeated));
assert(isnumeric(intersection_points));

% Check lengths
assert(isequal(size(min_cut),[1 1]));
assert(length(boundaryEngagedAtMinCut) == length(indices_repeated(:,1)));
assert(length(indices_repeated(1,:))==1);
assert(length(indices_repeated(:,1))>=1);
assert(length(intersection_points(1,:))==2);
assert(length(intersection_points(:,1))==length(vertices(:,1)));

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));

%% Compare speeds of pre-calculation versus post-calculation versus a fast variant
fig_num = 9903;
figure(fig_num);
close(fig_num);

Niterations = 100;

vertices = [0 0; 2 0; 1 2; 0 1; 0 0]*5;
[unit_normal_vectors, unit_vertex_projection_vectors, vector_direction_of_unit_cut, flag_vertexIsNonConvex]  = fcn_VSkel_polytopeFindUnitDirectionVectors(vertices,-1);
max_edge_cuts = fcn_VSkel_polytopeFindMaxEdgeCut(vertices, unit_normal_vectors, unit_vertex_projection_vectors, (-1));
[sphereRadii, definingBoundaries] = fcn_VSkel_polytopeFindEnclosedSpheres(vertices, unit_normal_vectors, unit_vertex_projection_vectors, vector_direction_of_unit_cut, flag_vertexIsNonConvex, max_edge_cuts, (-1));



% Do calculation without pre-calculation
tic;
for ith_test = 1:Niterations
    [min_cut, boundaryEngagedAtMinCut, indices_repeated, intersection_points] = ...
        fcn_VSkel_polytopeFindMinimumEnclosedSphere(vertices, ...
        sphereRadii, definingBoundaries, unit_normal_vectors, unit_vertex_projection_vectors, ...
        vector_direction_of_unit_cut, ([]));
end
slow_method = toc;

% Do calculation with pre-calculation, FAST_MODE on
tic;
for ith_test = 1:Niterations
    [min_cut, boundaryEngagedAtMinCut, indices_repeated, intersection_points] = ...
        fcn_VSkel_polytopeFindMinimumEnclosedSphere(vertices, ...
        sphereRadii, definingBoundaries, unit_normal_vectors, unit_vertex_projection_vectors, ...
        vector_direction_of_unit_cut, (-1));
end
fast_method = toc;

% Plot results as bar chart
figure(373737);
clf;
X = categorical({'Normal mode','Fast mode'});
X = reordercats(X,{'Normal mode','Fast mode'}); % Forces bars to appear in this exact order, not alphabetized
Y = [slow_method fast_method ]*1000/Niterations;
bar(X,Y)
ylabel('Execution time (Milliseconds)')


% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));
