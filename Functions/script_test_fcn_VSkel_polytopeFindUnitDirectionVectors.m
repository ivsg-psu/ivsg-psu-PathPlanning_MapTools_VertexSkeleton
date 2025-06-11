% script_test_fcn_VSkel_polytopeFindUnitDirectionVectors
% Tests function: fcn_VSkel_polytopeFindUnitDirectionVectors

% REVISION HISTORY:
% 2022_02_15
% -- first written by S. Brennan
% 2025_04_29 by Sean Brennan
% -- updated main script for VSkel library implementation
% 2025_05_03 by Sean Brennan
% -- updated main script to separate demonstration examples, basic
% examples, and fast tests
% 2025_05_14 by Sean Brennan
% -- added case where calcualtions can include line segments
% 2025_05_14 by Sean Brennan
% -- added case where calcualtions can include one point



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

%% Demonstration case 1: random 2D polytopes
fig_num = 0001;
figure(fig_num);
clf;

clear vertices
vertices{1} = [0 0; 1 0; 1 1]*5;
vertices{2} = [2 2; 3 4; 1 5]*5;
polytopeStructure = fcn_VSkel_polytopeFillStructureFromVertices(vertices, (-1));
Nvertices = length(polytopeStructure.polyPatch.Vertices(:,1));

% Call the function
[unit_normal_vectors, vector_direction_of_unit_cut, flag_vertexIsNonConvex]  = ...
    fcn_VSkel_polytopeFindUnitDirectionVectors(polytopeStructure,fig_num);

% Check variable types
assert(length(unit_normal_vectors(:,1)) == Nvertices);

assert(length(vector_direction_of_unit_cut(:,1)) == Nvertices);
assert(length(flag_vertexIsNonConvex(:,1)) == Nvertices);

% Check that all unit vectors are unit length
unit_normal_vectors_length = sum(unit_normal_vectors.^2,2).^0.5;
assert(all(abs(unit_normal_vectors_length - 1)<1E-10)==true);

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
vertices = [0 0; 2 0; 1 2; 0 1]*5;
polytopeStructure = fcn_VSkel_polytopeFillStructureFromVertices(vertices, (-1));
Nvertices = length(polytopeStructure.polyPatch.Vertices(:,1));

% Call the function
[unit_normal_vectors, vector_direction_of_unit_cut, flag_vertexIsNonConvex]  = ...
    fcn_VSkel_polytopeFindUnitDirectionVectors(polytopeStructure,fig_num);

% Check variable types
assert(length(unit_normal_vectors(:,1)) == Nvertices);

assert(length(vector_direction_of_unit_cut(:,1)) == Nvertices);
assert(length(flag_vertexIsNonConvex(:,1)) == Nvertices);

% Check that all unit vectors are unit length
unit_normal_vectors_length = sum(unit_normal_vectors.^2,2).^0.5;
assert(all(abs(unit_normal_vectors_length - 1)<1E-10)==true);

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Example case 2: square
fig_num = 1002;
figure(fig_num);
clf;

vertices = [0 0; 1 0; 1 1; 0 1]*10;
polytopeStructure = fcn_VSkel_polytopeFillStructureFromVertices(vertices, (-1));
Nvertices = length(polytopeStructure.polyPatch.Vertices(:,1));

% Call the function
[unit_normal_vectors, vector_direction_of_unit_cut, flag_vertexIsNonConvex]  = ...
    fcn_VSkel_polytopeFindUnitDirectionVectors(polytopeStructure,fig_num);

% Check variable types
assert(length(unit_normal_vectors(:,1)) == Nvertices);

assert(length(vector_direction_of_unit_cut(:,1)) == Nvertices);
assert(length(flag_vertexIsNonConvex(:,1)) == Nvertices);

% Check that all unit vectors are unit length
unit_normal_vectors_length = sum(unit_normal_vectors.^2,2).^0.5;
assert(all(abs(unit_normal_vectors_length - 1)<1E-10)==true);

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Example case 3: wide rectangle
fig_num = 1003;
figure(fig_num);
clf;

vertices = [0 0; 1 0; 1 0.5; 0 0.5]*10;
polytopeStructure = fcn_VSkel_polytopeFillStructureFromVertices(vertices, (-1));
Nvertices = length(polytopeStructure.polyPatch.Vertices(:,1));

% Call the function
[unit_normal_vectors, vector_direction_of_unit_cut, flag_vertexIsNonConvex]  = ...
    fcn_VSkel_polytopeFindUnitDirectionVectors(polytopeStructure,fig_num);

% Check variable types
assert(length(unit_normal_vectors(:,1)) == Nvertices);

assert(length(vector_direction_of_unit_cut(:,1)) == Nvertices);
assert(length(flag_vertexIsNonConvex(:,1)) == Nvertices);

% Check that all unit vectors are unit length
unit_normal_vectors_length = sum(unit_normal_vectors.^2,2).^0.5;
assert(all(abs(unit_normal_vectors_length - 1)<1E-10)==true);

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Example case 4: tall rectangle
fig_num = 1004;
figure(fig_num);
clf;

vertices = [0 0; 0.5 0; 0.5 1; 0 1]*10;
polytopeStructure = fcn_VSkel_polytopeFillStructureFromVertices(vertices, (-1));
Nvertices = length(polytopeStructure.polyPatch.Vertices(:,1));

% Call the function
[unit_normal_vectors, vector_direction_of_unit_cut, flag_vertexIsNonConvex]  = ...
    fcn_VSkel_polytopeFindUnitDirectionVectors(polytopeStructure,fig_num);

% Check variable types
assert(length(unit_normal_vectors(:,1)) == Nvertices);

assert(length(vector_direction_of_unit_cut(:,1)) == Nvertices);
assert(length(flag_vertexIsNonConvex(:,1)) == Nvertices);

% Check that all unit vectors are unit length
unit_normal_vectors_length = sum(unit_normal_vectors.^2,2).^0.5;
assert(all(abs(unit_normal_vectors_length - 1)<1E-10)==true);

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Example case 5: goofy polytope
fig_num = 1005;
figure(fig_num);
clf;

vertices = [0 0; 10 0; 5 15; 4 17; 1 13; 0 5];
polytopeStructure = fcn_VSkel_polytopeFillStructureFromVertices(vertices, (-1));
Nvertices = length(polytopeStructure.polyPatch.Vertices(:,1));

% Call the function
[unit_normal_vectors, vector_direction_of_unit_cut, flag_vertexIsNonConvex]  = ...
    fcn_VSkel_polytopeFindUnitDirectionVectors(polytopeStructure,fig_num);

% Check variable types
assert(length(unit_normal_vectors(:,1)) == Nvertices);

assert(length(vector_direction_of_unit_cut(:,1)) == Nvertices);
assert(length(flag_vertexIsNonConvex(:,1)) == Nvertices);

% Check that all unit vectors are unit length
unit_normal_vectors_length = sum(unit_normal_vectors.^2,2).^0.5;
assert(all(abs(unit_normal_vectors_length - 1)<1E-10)==true);

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Basic example of vertex calculation - non-convex 2D polytope
fig_num = 1006;
figure(fig_num);
clf;

% this polytope has a vertical wall
vertices = [0 0; 3 0; 5 5; 7  0; 10 0; 5 10; 0 5];
polytopeStructure = fcn_VSkel_polytopeFillStructureFromVertices(vertices, (-1));
Nvertices = length(polytopeStructure.polyPatch.Vertices(:,1));

% Call the function
[unit_normal_vectors, vector_direction_of_unit_cut, flag_vertexIsNonConvex]  = ...
    fcn_VSkel_polytopeFindUnitDirectionVectors(polytopeStructure,fig_num);

% Check variable types
assert(length(unit_normal_vectors(:,1)) == Nvertices);

assert(length(vector_direction_of_unit_cut(:,1)) == Nvertices);
assert(length(flag_vertexIsNonConvex(:,1)) == Nvertices);

% Check that all unit vectors are unit length
unit_normal_vectors_length = sum(unit_normal_vectors.^2,2).^0.5;
assert(all(abs(unit_normal_vectors_length - 1)<1E-10)==true);

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Basic example of vertex calculation - line segment
fig_num = 1007;
figure(fig_num);
clf;

% each face of this polytope has only 2 points - e.g. it is a line segment
clear vertices
vertices{1} = [0 0; 4 0];
vertices{2} = [1 1; 1 5];

polytopeStructure = fcn_VSkel_polytopeFillStructureFromVertices(vertices, (-1));
Nvertices = length(polytopeStructure.polyPatch.Vertices(:,1));

% Call the function
[unit_normal_vectors, vector_direction_of_unit_cut, flag_vertexIsNonConvex]  = ...
    fcn_VSkel_polytopeFindUnitDirectionVectors(polytopeStructure,fig_num);

% Check variable types
assert(length(unit_normal_vectors(:,1)) == Nvertices);

assert(length(vector_direction_of_unit_cut(:,1)) == Nvertices);
assert(length(flag_vertexIsNonConvex(:,1)) == Nvertices);

% Check that all unit vectors are unit length
unit_normal_vectors_length = sum(unit_normal_vectors.^2,2).^0.5;
assert(all(abs(unit_normal_vectors_length - 1)<1E-10)==true);

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));



%% Basic example of vertex calculation - single point
fig_num = 1008;
figure(fig_num);
clf;

% each face of this polytope has only 1 points - e.g. it is a degenerate line segment
vertices{1} = [4 0];
vertices{2} = [1 1];
polytopeStructure = fcn_VSkel_polytopeFillStructureFromVertices(vertices, (-1));
Nvertices = length(polytopeStructure.polyPatch.Vertices(:,1));

% Call the function
[unit_normal_vectors, vector_direction_of_unit_cut, flag_vertexIsNonConvex]  = ...
    fcn_VSkel_polytopeFindUnitDirectionVectors(polytopeStructure,fig_num);

% Check variable types
assert(length(unit_normal_vectors(:,1)) == Nvertices);

assert(length(vector_direction_of_unit_cut(:,1)) == Nvertices);
assert(length(flag_vertexIsNonConvex(:,1)) == Nvertices);

% Check that all unit vectors are unit length
unit_normal_vectors_length = sum(unit_normal_vectors.^2,2).^0.5;
assert(all(abs(unit_normal_vectors_length - 1)<1E-10)==false);

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));



%% Basic testing examples - Enclosing 2D Polytopes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  ____            _        _______        _   _                ______                           _
% |  _ \          (_)      |__   __|      | | (_)              |  ____|                         | |
% | |_) | __ _ ___ _  ___     | | ___  ___| |_ _ _ __   __ _   | |__  __  ____ _ _ __ ___  _ __ | | ___  ___
% |  _ < / _` / __| |/ __|    | |/ _ \/ __| __| | '_ \ / _` |  |  __| \ \/ / _` | '_ ` _ \| '_ \| |/ _ \/ __|
% | |_) | (_| \__ \ | (__     | |  __/\__ \ |_| | | | | (_| |  | |____ >  < (_| | | | | | | |_) | |  __/\__ \
% |____/ \__,_|___/_|\___|    |_|\___||___/\__|_|_| |_|\__, |  |______/_/\_\__,_|_| |_| |_| .__/|_|\___||___/
%                                                       __/ |                             | |
%                                                      |___/                              |_|
%      __     ______            _           _               ___  _____     _____      _       _
%      \ \   |  ____|          | |         (_)             |__ \|  __ \   |  __ \    | |     | |
%  _____\ \  | |__   _ __   ___| | ___  ___ _ _ __   __ _     ) | |  | |  | |__) |__ | |_   _| |_ ___  _ __   ___  ___
% |______> > |  __| | '_ \ / __| |/ _ \/ __| | '_ \ / _` |   / /| |  | |  |  ___/ _ \| | | | | __/ _ \| '_ \ / _ \/ __|
%       / /  | |____| | | | (__| | (_) \__ \ | | | | (_| |  / /_| |__| |  | |  | (_) | | |_| | || (_) | |_) |  __/\__ \
%      /_/   |______|_| |_|\___|_|\___/|___/_|_| |_|\__, | |____|_____/   |_|   \___/|_|\__, |\__\___/| .__/ \___||___/
%                                                    __/ |                               __/ |        | |
%                                                   |___/                               |___/         |_|
% See: http://patorjk.com/software/taag/#p=display&v=0&f=Big&t=Basic%20Testing%20%20Examples%20%0A-%3E%20Enclosing%202D%20%20Polytopes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Basic example of vertex calculation - non-convex 2D polytope
fig_num = 2001;
figure(fig_num);
clf;

% this polytope is a rectangle embedded within a rectangle
clear vertices;
vertices{1} = [-2 -2; 12 -2; 12 7; -2 7];
vertices{2} = flipud([2 2; 8 2; 8 3; 2 3]);
polytopeStructure = fcn_VSkel_polytopeFillStructureFromVertices(vertices, (-1));
Nvertices = length(polytopeStructure.polyPatch.Vertices(:,1));

% Call the function
[unit_normal_vectors, vector_direction_of_unit_cut, flag_vertexIsNonConvex]  = ...
    fcn_VSkel_polytopeFindUnitDirectionVectors(polytopeStructure,fig_num);

% Check variable types
assert(length(unit_normal_vectors(:,1)) == Nvertices);

assert(length(vector_direction_of_unit_cut(:,1)) == Nvertices);
assert(length(flag_vertexIsNonConvex(:,1)) == Nvertices);

% Check that all unit vectors are unit length
unit_normal_vectors_length = sum(unit_normal_vectors.^2,2).^0.5;
assert(all(abs(unit_normal_vectors_length - 1)<1E-10)==true);

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));



%% Basic testing examples in 3D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  ____            _        _______        _   _                ______                           _                       ____  _____
% |  _ \          (_)      |__   __|      | | (_)              |  ____|                         | |                     |___ \|  __ \
% | |_) | __ _ ___ _  ___     | | ___  ___| |_ _ _ __   __ _   | |__  __  ____ _ _ __ ___  _ __ | | ___  ___    ______    __) | |  | |
% |  _ < / _` / __| |/ __|    | |/ _ \/ __| __| | '_ \ / _` |  |  __| \ \/ / _` | '_ ` _ \| '_ \| |/ _ \/ __|  |______|  |__ <| |  | |
% | |_) | (_| \__ \ | (__     | |  __/\__ \ |_| | | | | (_| |  | |____ >  < (_| | | | | | | |_) | |  __/\__ \            ___) | |__| |
% |____/ \__,_|___/_|\___|    |_|\___||___/\__|_|_| |_|\__, |  |______/_/\_\__,_|_| |_| |_| .__/|_|\___||___/           |____/|_____/
%                                                       __/ |                             | |
%                                                      |___/                              |_|
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Basic%20Testing%20%20Examples%20%20-%203D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Demonstration case: 3D polytope (cube)
fig_num = 3001;
figure(fig_num);
close(fig_num);

% The following is motivated by a square cube. It has 6
% external faces, and each vertex connects 3 faces.
clear vertices
vertices{1} = flipud([0 0 0; 0 1 0; 1 1 0; 1 0 0])*5;  % Bottom XY face
vertices{2} = flipud([0 0 0; 0 0 1; 0 1 1; 0 1 0])*5;  % Bottom YZ face
vertices{3} = flipud([0 0 0; 1 0 0; 1 0 1; 0 0 1])*5;  % Bottom XZ face
vertices{4} = flipud([1 1 1; 0 1 1; 0 0 1; 1 0 1])*5;  % Top XY face
vertices{5} = flipud([1 1 1; 1 0 1; 1 0 0; 1 1 0])*5;  % Top YZ face
vertices{6} = flipud([1 1 1; 1 1 0; 0 1 0; 0 1 1])*5;  % Top XZ face

polytopeStructure = fcn_VSkel_polytopeFillStructureFromVertices(vertices, -1);
Nvertices = length(polytopeStructure.polyPatch.Vertices(:,1));

% Call the function
[unit_normal_vectors, vector_direction_of_unit_cut, flag_vertexIsNonConvex]  = ...
    fcn_VSkel_polytopeFindUnitDirectionVectors(polytopeStructure,fig_num);

% Check variable types
Nfaces    = length(vertices);
assert(length(unit_normal_vectors(:,1)) == Nfaces);

assert(length(vector_direction_of_unit_cut(:,1)) == Nvertices);
assert(length(flag_vertexIsNonConvex(:,1)) == Nvertices);

% Check that all unit vectors are unit length
unit_normal_vectors_length = sum(unit_normal_vectors.^2,2).^0.5;
assert(all(abs(unit_normal_vectors_length - 1)<1E-10)==true);

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));


%% Demonstration case: 3D polytope (square bottom pyramid)
%                    V5
%                   /=\\
%                  /===\ \
%                 /=====\' \
%                /=======\'' \
%               /=========\ ' '\
%              /===========\''   \
%             /=============\ ' '  \
%            /===============\  F3'  \
%           /=======F2 =======\' ' ' ' \
%          /===================\' ' '  ' \
%         /=====================\' '   ' ' V3
%        /=======================\  '   ' /
%       /=========================\   ' /
%      /===========================\'  /
%     V1============================V2
%    

fig_num = 3002;
figure(fig_num);
close(fig_num);

clear vertices
vertices{1} = flipud([0 0 0; 0 2 0; 2 2 0; 2 0 0])*5;  % Bottom XY face
vertices{2} = flipud([0 0 0; 2 0 0; 1 1 2])*5;         % front face
vertices{3} = flipud([2 0 0; 2 2 0; 1 1 2])*5;         % right face
vertices{4} = flipud([2 2 0; 0 2 0; 1 1 2])*5;         % back face
vertices{5} = flipud([0 2 0; 0 0 0; 1 1 2])*5;         % left face

polytopeStructure = fcn_VSkel_polytopeFillStructureFromVertices(vertices, -1);
Nvertices = length(polytopeStructure.polyPatch.Vertices(:,1));

% Call the function
[unit_normal_vectors, vector_direction_of_unit_cut, flag_vertexIsNonConvex]  = ...
    fcn_VSkel_polytopeFindUnitDirectionVectors(polytopeStructure,fig_num);

% Check variable types
Nfaces    = length(vertices);
assert(length(unit_normal_vectors(:,1)) == Nfaces);

assert(length(vector_direction_of_unit_cut(:,1)) == Nvertices);
assert(length(flag_vertexIsNonConvex(:,1)) == Nvertices);

% Check that all unit vectors are unit length
unit_normal_vectors_length = sum(unit_normal_vectors.^2,2).^0.5;
assert(all(abs(unit_normal_vectors_length - 1)<1E-10)==true);

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Demonstration case: 3D polytope (triangular base pyramid)
%                  4
%                 / \ >
%                /   \  >
%               / (3) \   > 
%              /    1  \ 2 < > 3 
%             /      <  \   /   
%            /   <  (4)  \ /   
%           1-<-----------2
%
%
% Correct answer should be:
% 1     2     3     4
     
fig_num = 3003;
figure(fig_num);
close(fig_num);

clear vertices
vertices{1} = flipud([0 0 0; 1 1 0; 2 0 0])*5;    % Bottom XY face
vertices{2} = flipud([0 0 0; 2 0 0; 1 0.5 2])*5;  % front face
vertices{3} = flipud([2 0 0; 1 1 0; 1 0.5 2])*5;  % right face
vertices{4} = flipud([1 1 0; 0 0 0; 1 0.5 2])*5;  % left face


polytopeStructure = fcn_VSkel_polytopeFillStructureFromVertices(vertices, -1);
Nvertices = length(polytopeStructure.polyPatch.Vertices(:,1));

% Call the function
[unit_normal_vectors, vector_direction_of_unit_cut, flag_vertexIsNonConvex]  = ...
    fcn_VSkel_polytopeFindUnitDirectionVectors(polytopeStructure,fig_num);

% Check variable types
Nfaces    = length(vertices);
assert(length(unit_normal_vectors(:,1)) == Nfaces);

assert(length(vector_direction_of_unit_cut(:,1)) == Nvertices);
assert(length(flag_vertexIsNonConvex(:,1)) == Nvertices);

% Check that all unit vectors are unit length
unit_normal_vectors_length = sum(unit_normal_vectors.^2,2).^0.5;
assert(all(abs(unit_normal_vectors_length - 1)<1E-10)==true);

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Demonstration case: 3D polytope (cube inside out)
fig_num = 3004;
figure(fig_num);
close(fig_num);

% The following is motivated by a square cube. It has 6
% external faces, and each vertex connects 3 faces.
clear vertices
vertices{1} = ([0 0 0; 0 1 0; 1 1 0; 1 0 0])*5;  % Bottom XY face
vertices{2} = ([0 0 0; 0 0 1; 0 1 1; 0 1 0])*5;  % Bottom YZ face
vertices{3} = ([0 0 0; 1 0 0; 1 0 1; 0 0 1])*5;  % Bottom XZ face
vertices{4} = ([1 1 1; 0 1 1; 0 0 1; 1 0 1])*5;  % Top XY face
vertices{5} = ([1 1 1; 1 0 1; 1 0 0; 1 1 0])*5;  % Top YZ face
vertices{6} = ([1 1 1; 1 1 0; 0 1 0; 0 1 1])*5;  % Top XZ face

polytopeStructure = fcn_VSkel_polytopeFillStructureFromVertices(vertices, -1);
Nvertices = length(polytopeStructure.polyPatch.Vertices(:,1));

% Call the function
[unit_normal_vectors, vector_direction_of_unit_cut, flag_vertexIsNonConvex]  = ...
    fcn_VSkel_polytopeFindUnitDirectionVectors(polytopeStructure,fig_num);

% Check variable types
Nfaces    = length(vertices);
assert(length(unit_normal_vectors(:,1)) == Nfaces);

assert(length(vector_direction_of_unit_cut(:,1)) == Nvertices);
assert(length(flag_vertexIsNonConvex(:,1)) == Nvertices);

% Check that all unit vectors are unit length
unit_normal_vectors_length = sum(unit_normal_vectors.^2,2).^0.5;
assert(all(abs(unit_normal_vectors_length - 1)<1E-10)==true);

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Demonstration case: 3D polytope (cube inside cube)
fig_num = 3005;
figure(fig_num);
close(fig_num);

% The following is motivated by a square cube. It has 6
% external faces, and each vertex connects 3 faces.
clear vertices
vertices{1} = ([0 0 0; 0 1 0; 1 1 0; 1 0 0])*5 + ones(4,1)*[2.5 2.5 2.5];  % Bottom XY face
vertices{2} = ([0 0 0; 0 0 1; 0 1 1; 0 1 0])*5 + ones(4,1)*[2.5 2.5 2.5];  % Bottom YZ face
vertices{3} = ([0 0 0; 1 0 0; 1 0 1; 0 0 1])*5 + ones(4,1)*[2.5 2.5 2.5];  % Bottom XZ face
vertices{4} = ([1 1 1; 0 1 1; 0 0 1; 1 0 1])*5 + ones(4,1)*[2.5 2.5 2.5];  % Top XY face
vertices{5} = ([1 1 1; 1 0 1; 1 0 0; 1 1 0])*5 + ones(4,1)*[2.5 2.5 2.5];  % Top YZ face
vertices{6} = ([1 1 1; 1 1 0; 0 1 0; 0 1 1])*5 + ones(4,1)*[2.5 2.5 2.5];  % Top XZ face
vertices{7} = flipud([0 0 0; 0 1 0; 1 1 0; 1 0 0])*10;  % Bottom XY face
vertices{8} = flipud([0 0 0; 0 0 1; 0 1 1; 0 1 0])*10;  % Bottom YZ face
vertices{9} = flipud([0 0 0; 1 0 0; 1 0 1; 0 0 1])*10;  % Bottom XZ face
vertices{10} = flipud([1 1 1; 0 1 1; 0 0 1; 1 0 1])*10;  % Top XY face
vertices{11} = flipud([1 1 1; 1 0 1; 1 0 0; 1 1 0])*10;  % Top YZ face
vertices{12} = flipud([1 1 1; 1 1 0; 0 1 0; 0 1 1])*10;  % Top XZ face

polytopeStructure = fcn_VSkel_polytopeFillStructureFromVertices(vertices, -1);
Nvertices = length(polytopeStructure.polyPatch.Vertices(:,1));

% Call the function
[unit_normal_vectors, vector_direction_of_unit_cut, flag_vertexIsNonConvex]  = ...
    fcn_VSkel_polytopeFindUnitDirectionVectors(polytopeStructure,fig_num);

% Check variable types
Nfaces    = length(vertices);
assert(length(unit_normal_vectors(:,1)) == Nfaces);

assert(length(vector_direction_of_unit_cut(:,1)) == Nvertices);
assert(length(flag_vertexIsNonConvex(:,1)) == Nvertices);

% Check that all unit vectors are unit length
unit_normal_vectors_length = sum(unit_normal_vectors.^2,2).^0.5;
assert(all(abs(unit_normal_vectors_length - 1)<1E-10)==true);

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));


%% Demonstration case: 3D polytope (cube created by adjacent faces)
fig_num = 3006;
figure(fig_num);
close(fig_num);

% The following is motivated by a square cube. It has 6
% external faces, and each vertex connects 3 faces.
clear vertices
vertices{1} = flipud([0 0 0; 0 1 0; 1 1 0])*5;  % Bottom XY face (1 of 2)
vertices{2} = flipud([0 0 0; 1 1 0; 1 0 0])*5;  % Bottom XY face (2 of 2)
vertices{3} = flipud([0 0 0; 0 0 1; 0 1 1; 0 1 0])*5;  % Bottom YZ face
vertices{4} = flipud([0 0 0; 1 0 0; 1 0 1; 0 0 1])*5;  % Bottom XZ face
vertices{5} = flipud([1 1 1; 0 1 1; 0 0 1; 1 0 1])*5;  % Top XY face
vertices{6} = flipud([1 1 1; 1 0 1; 1 0 0; 1 1 0])*5;  % Top YZ face
vertices{7} = flipud([1 1 1; 1 1 0; 0 1 0; 0 1 1])*5;  % Top XZ face

polytopeStructure = fcn_VSkel_polytopeFillStructureFromVertices(vertices, -1);
Nvertices = length(polytopeStructure.polyPatch.Vertices(:,1));

% Call the function
[unit_normal_vectors, vector_direction_of_unit_cut, flag_vertexIsNonConvex]  = ...
    fcn_VSkel_polytopeFindUnitDirectionVectors(polytopeStructure,fig_num);

% Check variable types
Nfaces    = length(vertices);
assert(length(unit_normal_vectors(:,1)) == Nfaces);

assert(length(vector_direction_of_unit_cut(:,1)) == Nvertices);
assert(length(flag_vertexIsNonConvex(:,1)) == Nvertices);

% Check that all unit vectors are unit length
unit_normal_vectors_length = sum(unit_normal_vectors.^2,2).^0.5;
assert(all(abs(unit_normal_vectors_length - 1)<1E-10)==true);

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Demonstration case: 3D polytope (shape created by adjoint patches)
fig_num = 3007;
figure(fig_num);
close(fig_num);

% The following is motivated by a square cube. It has 6
% external faces, and each vertex connects 3 faces.
clear vertices
vertices = cell(6,1);
vertices{1} = ([0 0 0; 10 0 0; 10 10 0; 0 10 0]);  % Bottom XY face
vertices{2} = ([10 0 0; 7.5 0 2.5; 7.5 10 2.5; 10 10 0]);  % Right YZ face on base of bottom 
vertices{3} = ([2.5 0 2.5; 2.5 10 2.5; 7.5 10 2.5; 7.5 0 2.5]);  % Top XY face on base of bottom
vertices{4} = ([2.5 0 2.5; 0 0 0; 0 10 0; 2.5 10 2.5]);  % Left YZ face on inset of bottom
vertices{5} = ([0 0 0; 2.5 0 2.5; 7.5 0 2.5; 10 0 0]);  % Front XZ of bottom
vertices{6} = ([0 10 0; 10 10 0; 7.5 10 2.5; 2.5 10 2.5]);  % Front XZ of bottom

% Generate the next side by rotation
% Define the angle of rotation in degrees
angle = 90; % Rotate by 45 degrees

% Convert angle to radians
theta = deg2rad(angle);

% Create the rotation matrix for rotation about the y-axis
R = [cos(theta) 0 sin(theta);
     0         1 0;
     -sin(theta) 0 cos(theta)];

rotatedPoints = cell(length(vertices),1);
for ith_face = 1:length(vertices)
    points = vertices{ith_face};
    Npoints = length(points(:,1));
    points_shifted = points - ones(Npoints,1)*[5 0 5];

    % Rotate the points
    rotated_shifted_points = (R * points_shifted')'; % Transpose for matrix multiplication

    rotatedPoints{ith_face} = rotated_shifted_points + ones(Npoints,1)*[5 0 5];
end

full_vertices = [vertices; rotatedPoints];

polytopeStructure = fcn_VSkel_polytopeFillStructureFromVertices(full_vertices, -1);
Nvertices = length(polytopeStructure.polyPatch.Vertices(:,1));

% Call the function
[unit_normal_vectors, vector_direction_of_unit_cut, flag_vertexIsNonConvex]  = ...
    fcn_VSkel_polytopeFindUnitDirectionVectors(polytopeStructure,fig_num);

% Check variable types
Nfaces    = length(full_vertices);
assert(length(unit_normal_vectors(:,1)) == Nfaces);

assert(length(vector_direction_of_unit_cut(:,1)) == Nvertices);
assert(length(flag_vertexIsNonConvex(:,1)) == Nvertices);

% Check that all unit vectors are unit length
unit_normal_vectors_length = sum(unit_normal_vectors.^2,2).^0.5;
assert(all(abs(unit_normal_vectors_length - 1)<1E-10)==true);

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Demonstration case: 3D polytope (triangular base pyramid with cut-out on face 1 and 2)
%                  4
%                 / \ >
%                /   \  >
%               / (3) \   > 
%              /    1  \ 2 < > 3 
%             /      <  \   /   
%            /   <  (4)  \ /   
%           1-<-----------2
%
%
% Correct answer should be:
% 1     2     3     4
     
fig_num = 3008;
figure(fig_num);
close(fig_num);

clear vertices
vertices{1} = flipud([0 0 0; 1 2 0; 2 0 0; 1 1 0])*5;    % Bottom XY face with cut-out
vertices{2} = flipud([0 0 0; 1 1 0; 1 0.5 2])*5;  % front face left cut out
vertices{3} = flipud([1 1 0; 2 0 0; 1 0.5 2])*5;  % front face right cut out
vertices{4} = flipud([2 0 0; 1 2 0; 1 0.5 2])*5;  % right face
vertices{5} = flipud([1 2 0; 0 0 0; 1 0.5 2])*5;  % left face


polytopeStructure = fcn_VSkel_polytopeFillStructureFromVertices(vertices, -1);
Nvertices = length(polytopeStructure.polyPatch.Vertices(:,1));

% Call the function
[unit_normal_vectors, vector_direction_of_unit_cut, flag_vertexIsNonConvex]  = ...
    fcn_VSkel_polytopeFindUnitDirectionVectors(polytopeStructure,fig_num);

% Check variable types
Nfaces    = length(vertices);
assert(length(unit_normal_vectors(:,1)) == Nfaces);

assert(length(vector_direction_of_unit_cut(:,1)) == Nvertices);
assert(length(flag_vertexIsNonConvex(:,1)) == Nvertices);

% Check that all unit vectors are unit length
unit_normal_vectors_length = sum(unit_normal_vectors.^2,2).^0.5;
assert(all(abs(unit_normal_vectors_length - 1)<1E-10)==true);

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Demonstration case: 3D polytope (square bottom pyramid cut into cube)
%                    V5
%                   /=\\
%                  /===\ \
%                 /=====\' \
%                /=======\'' \
%               /=========\ ' '\
%              /===========\''   \
%             /=============\ ' '  \
%            /===============\  F3'  \
%           /=======F2 =======\' ' ' ' \
%          /===================\' ' '  ' \
%         /=====================\' '   ' ' V3
%        /=======================\  '   ' /
%       /=========================\   ' /
%      /===========================\'  /
%     V1============================V2
%    

URHERE
fig_num = 3009;
figure(fig_num);
close(fig_num);

clear vertices
% vertices{1} = flipud([0 0 0; 0 2 0; 2 2 0; 2 0 0])*5;  % Pyramid bottom XY face
vertices{1} = ([0 0 0; 2 0 0; 1 1 2])*5;         % Pyramid front face
vertices{2} = ([2 0 0; 2 2 0; 1 1 2])*5;         % Pyramid right face
vertices{3} = ([2 2 0; 0 2 0; 1 1 2])*5;         % Pyramid back face
vertices{4} = ([0 2 0; 0 0 0; 1 1 2])*5;         % Pyramid left face

vertices{5} = flipud([0 0 0; 0 0 1; 0 1 1; 0 1 0])*10;  % Cube bottom YZ face
vertices{6} = flipud([0 0 0; 1 0 0; 1 0 1; 0 0 1])*10;  % Cube bottom XZ face
vertices{7} = flipud([1 1 1; 0 1 1; 0 0 1; 1 0 1])*10;  % Cube top XY face
vertices{8} = flipud([1 1 1; 1 0 1; 1 0 0; 1 1 0])*10;  % Cube top YZ face
vertices{9} = flipud([1 1 1; 1 1 0; 0 1 0; 0 1 1])*10;  % Cube top XZ face

polytopeStructure = fcn_VSkel_polytopeFillStructureFromVertices(vertices, -1);
Nvertices = length(polytopeStructure.polyPatch.Vertices(:,1));

% Call the function
[unit_normal_vectors, vector_direction_of_unit_cut, flag_vertexIsNonConvex]  = ...
    fcn_VSkel_polytopeFindUnitDirectionVectors(polytopeStructure,fig_num);

% Check variable types
Nfaces    = length(vertices);
assert(length(unit_normal_vectors(:,1)) == Nfaces);

assert(length(vector_direction_of_unit_cut(:,1)) == Nvertices);
assert(length(flag_vertexIsNonConvex(:,1)) == Nvertices);

% Check that all unit vectors are unit length
unit_normal_vectors_length = sum(unit_normal_vectors.^2,2).^0.5;
assert(all(abs(unit_normal_vectors_length - 1)<1E-10)==true);

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

vertices = [0 0; 2 0; 1 2; 0 1]*5;
polytopeStructure = fcn_VSkel_polytopeFillStructureFromVertices(vertices, (-1));
Nvertices = length(polytopeStructure.polyPatch.Vertices(:,1));

% Call the function
[unit_normal_vectors, vector_direction_of_unit_cut, flag_vertexIsNonConvex]  = ...
    fcn_VSkel_polytopeFindUnitDirectionVectors(polytopeStructure,[]);

% Check variable types
assert(length(unit_normal_vectors(:,1)) == Nvertices);

assert(length(vector_direction_of_unit_cut(:,1)) == Nvertices);
assert(length(flag_vertexIsNonConvex(:,1)) == Nvertices);

% Check that all unit vectors are unit length
unit_normal_vectors_length = sum(unit_normal_vectors.^2,2).^0.5;
assert(all(abs(unit_normal_vectors_length - 1)<1E-10)==true);

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));

%% Basic example of vertex calculation - non-normal wall shrinking, NO FIGURE, FAST MODE
fig_num = 9902;
figure(fig_num);
close(fig_num);

vertices = [0 0; 2 0; 1 2; 0 1]*5;
polytopeStructure = fcn_VSkel_polytopeFillStructureFromVertices(vertices, (-1));
Nvertices = length(polytopeStructure.polyPatch.Vertices(:,1));

% Call the function
[unit_normal_vectors, vector_direction_of_unit_cut, flag_vertexIsNonConvex]  = ...
    fcn_VSkel_polytopeFindUnitDirectionVectors(polytopeStructure,-1);

% Check variable types
assert(length(unit_normal_vectors(:,1)) == Nvertices);

assert(length(vector_direction_of_unit_cut(:,1)) == Nvertices);
assert(length(flag_vertexIsNonConvex(:,1)) == Nvertices);

% Check that all unit vectors are unit length
unit_normal_vectors_length = sum(unit_normal_vectors.^2,2).^0.5;
assert(all(abs(unit_normal_vectors_length - 1)<1E-10)==true);

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));

%% Compare speeds of pre-calculation versus post-calculation versus a fast variant
fig_num = 9903;
figure(fig_num);
close(fig_num);

Niterations = 100;

% Set up polytopes
% load('testData_fcn_VSkel_plotPolytopes.mat','polytopes','polytopes2');

vertices = [0 0; 2 0; 1 2; 0 1]*5;
polytopeStructure = fcn_VSkel_polytopeFillStructureFromVertices(vertices, (-1));

% Do calculation without pre-calculation
tic;
for ith_test = 1:Niterations
    [unit_normal_vectors, vector_direction_of_unit_cut, flag_vertexIsNonConvex]  = ...
        fcn_VSkel_polytopeFindUnitDirectionVectors(polytopeStructure,[]);
end
slow_method = toc;

% Do calculation with pre-calculation, FAST_MODE on
tic;
for ith_test = 1:Niterations
    [unit_normal_vectors, vector_direction_of_unit_cut, flag_vertexIsNonConvex]  = ...
        fcn_VSkel_polytopeFindUnitDirectionVectors(polytopeStructure,-1);
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
