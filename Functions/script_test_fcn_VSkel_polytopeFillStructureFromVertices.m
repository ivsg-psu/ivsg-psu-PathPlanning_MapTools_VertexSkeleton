% script_test_fcn_VSkel_polytopeFillStructureFromVertices
% Tests function: fcn_VSkel_polytopeFillStructureFromVertices

% REVISION HISTORY:
% 2025_05_24 by Sean Brennan
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

%% Demonstration case 1: multiple 2D polytope
fig_num = 0001;
figure(fig_num);
clf;

% Call the function
clear vertices
vertices{1} = [0 0; 1 0; 1 1];
vertices{2} = [2 2; 3 4; 1 5];
polytopeStructure = fcn_VSkel_polytopeFillStructureFromVertices(vertices, (fig_num));

% Check variable types
assert(isstruct(polytopeStructure));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));


%% Basic testing examples in 2D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  ____            _        _______        _   _                ______                           _                       ___  _____
% |  _ \          (_)      |__   __|      | | (_)              |  ____|                         | |                     |__ \|  __ \
% | |_) | __ _ ___ _  ___     | | ___  ___| |_ _ _ __   __ _   | |__  __  ____ _ _ __ ___  _ __ | | ___  ___    ______     ) | |  | |
% |  _ < / _` / __| |/ __|    | |/ _ \/ __| __| | '_ \ / _` |  |  __| \ \/ / _` | '_ ` _ \| '_ \| |/ _ \/ __|  |______|   / /| |  | |
% | |_) | (_| \__ \ | (__     | |  __/\__ \ |_| | | | | (_| |  | |____ >  < (_| | | | | | | |_) | |  __/\__ \            / /_| |__| |
% |____/ \__,_|___/_|\___|    |_|\___||___/\__|_|_| |_|\__, |  |______/_/\_\__,_|_| |_| |_| .__/|_|\___||___/           |____|_____/
%                                                       __/ |                             | |
%                                                      |___/                              |_|
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Basic%20Testing%20%20Examples%20%20-%202D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Basic case 1: single 2D polytope with 3 vertices
fig_num = 1001;
figure(fig_num);
clf;

% Call the function
vertices = [0 0; 1 0; 1 1];
polytopeStructure = fcn_VSkel_polytopeFillStructureFromVertices(vertices, (fig_num));

% Check variable types
assert(isstruct(polytopeStructure));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));



%% Basic case 2: multiple 2D polytope
fig_num = 1002;
figure(fig_num);
clf;

% Call the function
clear vertices
vertices{1} = [0 0; 1 0; 1 1];
vertices{2} = [2 2; 3 4; 1 5];
polytopeStructure = fcn_VSkel_polytopeFillStructureFromVertices(vertices, (fig_num));

% Check variable types
assert(isstruct(polytopeStructure));

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

% % Works in 3D also
% clear poly
% poly.Vertices = [0 0 0; 1 0 1; 1 1 0; 2 2 0; 3 2 0; 4 3 5];
% poly.Faces = [1 2; 2 3; 3 1; 4 5; 5 6; 6 4];
% poly.FaceVertexCData = [0; 1; 0.5];
% poly.FaceColor = 'none'; % flat
% poly.EdgeColor = 'flat';
% vertexNumbering = (1:length(poly.Vertices(:,1)))';
% colorIndices = mod(vertexNumbering-1,length(colorOrdering(:,1)))+1;
% poly.FaceVertexCData = colorOrdering(colorIndices,:);
% poly.LineWidth = 2;
% patch(poly);

fig_num = 3001;
figure(fig_num);
clf;

% Call the function
clear vertices
vertices{1} = [0 0 0; 1 0 1; 1 1 0];
vertices{2} = [2 2 0; 3 2 0; 4 3 5];
polytopeStructure = fcn_VSkel_polytopeFillStructureFromVertices(vertices, (fig_num));

% Check variable types
assert(isstruct(polytopeStructure));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));




%% Demonstration case: 3D polytope (square)
fig_num = 2002;
figure(fig_num);
close(fig_num);

% The following is motivated by a square cube. It has 6
% external faces, and each vertex connects 3 faces.
clear vertices
vertices{1} = [0 0 0; 1 0 0; 1 1 0; 0 1 0];
vertices{2} = [0 0 0; 0 1 0; 0 1 1; 0 0 1];
vertices{3} = [0 0 0; 0 0 1; 1 0 1; 1 0 0];
vertices{4} = [1 1 1; 0 1 1; 1 1 1; 1 0 1];
vertices{5} = [1 1 1; 1 0 1; 1 0 0; 1 1 0];
vertices{6} = [1 1 1; 1 1 0; 0 1 0; 0 1 1];


% Call the function
polytopeStructure = fcn_VSkel_polytopeFillStructureFromVertices(vertices, (fig_num));

% Check variable types
assert(isstruct(polytopeStructure));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));




%% Demonstration case: 3D polytope (square bottom pyramid)
%                    V1
%                   /=\\
%                  /===\ \
%                 /=====\' \
%                /=======\'' \
%               /=========\ ' '\
%              /===========\''   \
%             /=============\ ' '  \
%            /===============\  F2'  \
%           /=======F1 =======\' ' ' ' \
%          /===================\' ' '  ' \
%         /=====================\' '   ' ' V3
%        /=======================\  '   ' /
%       /=========================\   ' /
%      /===========================\'  /
%     V5============================V2
%
% Correct answer should be:
% 1     2     3     5
% 1     2     4     5
% 1     3     4     5
% 2     3     4     5
URHERE     

fig_num = 2002;
figure(fig_num);
close(fig_num);

% The following is motivated by a diamond shape in 3D, e.g. a square
% pyramid stacked atop another square pyramid. It has 6 vertices, 4 in the
% middle in a square pattern, and one on each side of the square. It has 8
% external faces, and each vertex connects 4 faces.
clear vertices
Nvertices = 5;
cell_array_edges_in_vertices{1,1} = [1 2 3 4];
cell_array_edges_in_vertices{2,1} = [1 2 5];
cell_array_edges_in_vertices{3,1} = [2 3 5];
cell_array_edges_in_vertices{4,1} = [3 4 5];
cell_array_edges_in_vertices{5,1} = [4 1 5];

Nfaces = 5;

% Call the function
polytopeStructure = fcn_VSkel_polytopeFillStructureFromVertices(vertices, (fig_num));

% Check variable types
assert(isstruct(polytopeStructure));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));



%% Demonstration case: 3D polytope (square cube)
%             3-------4
%            /  1    / |
%           /       /  |
%          2-------1 3 |8
%          |       |  /
%          |   2   | /
%          |       |/
%          6-------5   (4 on left, 5 on back, 6 on bottom)
%
% Correct answer should be:
% 1     2     3     4
% 1     2     3     5
% 1     2     3     6
% 1     2     4     5
% 1     2     4     6
% 1     3     4     5
% 1     3     5     6
% 1     4     5     6
% 2     3     4     6
% 2     3     5     6
% 2     4     5     6
% 3     4     5     6
     
fig_num = 2003;
figure(fig_num);
close(fig_num);

% The following is motivated by a diamond shape in 3D, e.g. a square
% pyramid stacked atop another square pyramid. It has 6 vertices, 4 in the
% middle in a square pattern, and one on each side of the square. It has 8
% external faces, and each vertex connects 4 faces.
clear vertices
Nvertices = 6;
cell_array_edges_in_vertices{1,1} = [1 2 3];
cell_array_edges_in_vertices{2,1} = [1 2 4];
cell_array_edges_in_vertices{3,1} = [1 4 5];
cell_array_edges_in_vertices{4,1} = [1 3 5];
cell_array_edges_in_vertices{5,1} = [2 3 6];
cell_array_edges_in_vertices{6,1} = [2 4 6];
cell_array_edges_in_vertices{7,1} = [4 5 6];
cell_array_edges_in_vertices{8,1} = [3 5 6];

Nfaces = 6;

% Call the function
polytopeStructure = fcn_VSkel_polytopeFillStructureFromVertices(vertices, (fig_num));

% Check variable types
assert(isstruct(polytopeStructure));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));



%% Demonstration case: 3D polytope (triangular base pyramid)
%                  1
%                 / \ >
%                /   \  >
%               / (3) \   > 
%              /    1  \ 2 < > 4 
%             /      <  \   /   
%            /   <  (4)  \ /   
%           2-<-----------3
%
%
% Correct answer should be:
% 1     2     3     4
     
fig_num = 2003;
figure(fig_num);
close(fig_num);

% The following is motivated by a diamond shape in 3D, e.g. a square
% pyramid stacked atop another square pyramid. It has 6 vertices, 4 in the
% middle in a square pattern, and one on each side of the square. It has 8
% external faces, and each vertex connects 4 faces.
clear vertices
Nvertices = 4;
cell_array_edges_in_vertices{1,1} = [1 2 3];
cell_array_edges_in_vertices{2,1} = [1 3 4];
cell_array_edges_in_vertices{3,1} = [1 2 4];
cell_array_edges_in_vertices{4,1} = [2 3 4];

Nfaces = 4;

% Call the function
polytopeStructure = fcn_VSkel_polytopeFillStructureFromVertices(vertices, (fig_num));

% Check variable types
assert(isstruct(polytopeStructure));

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

% Fill some vertices with edges. Each vertex is created by its own edge and
% the one prior to it
clear vertices
Nvertices = 30;
cell_array_edges_in_vertices = cell(Nvertices,1);
for ith_vertex = 1:Nvertices
    raw_edges = [ith_vertex-1; ith_vertex];
    cell_array_edges_in_vertices{ith_vertex,1} = mod(raw_edges,Nvertices)+1;
end

NE = Nvertices;

% Call the function
edge_permutations = fcn_VSkel_polytopeFillStructureFromVertices(cell_array_edges_in_vertices, NE, ([]));

% Check variable types
assert(isstruct(polytopeStructure));


% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));

%% Basic example of vertex calculation - non-normal wall shrinking, NO FIGURE, FAST MODE
fig_num = 9902;
figure(fig_num);
close(fig_num);


% Fill some vertices with edges. Each vertex is created by its own edge and
% the one prior to it
clear vertices
Nvertices = 30;
cell_array_edges_in_vertices = cell(Nvertices,1);
for ith_vertex = 1:Nvertices
    raw_edges = [ith_vertex-1; ith_vertex];
    cell_array_edges_in_vertices{ith_vertex,1} = mod(raw_edges,Nvertices)+1;
end

NE = Nvertices;

% Call the function
edge_permutations = fcn_VSkel_polytopeFillStructureFromVertices(cell_array_edges_in_vertices, NE, (-1));

% Check variable types
assert(isstruct(polytopeStructure));

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));

%% Compare speeds of pre-calculation versus post-calculation versus a fast variant
fig_num = 9903;
figure(fig_num);
close(fig_num);

Niterations = 100;


% Fill some vertices with edges. Each vertex is created by its own edge and
% the one prior to it
clear vertices
Nvertices = 30;
cell_array_edges_in_vertices = cell(Nvertices,1);
for ith_vertex = 1:Nvertices
    raw_edges = [ith_vertex-1; ith_vertex];
    cell_array_edges_in_vertices{ith_vertex,1} = mod(raw_edges,Nvertices)+1;
end

NE = Nvertices;

% Do calculation without pre-calculation
tic;
for ith_test = 1:Niterations
    % Call the function
    edge_permutations = fcn_VSkel_polytopeFillStructureFromVertices(cell_array_edges_in_vertices, NE, ([]));
end
slow_method = toc;

% Do calculation with pre-calculation, FAST_MODE on
tic;
for ith_test = 1:Niterations
    % Call the function
    edge_permutations = fcn_VSkel_polytopeFillStructureFromVertices(cell_array_edges_in_vertices, NE, (-1));
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
