% script_test_fcn_VSkel_findEdgePermutations
% Tests function: fcn_VSkel_findEdgePermutations

% REVISION HISTORY:
% 2025_05_18 by Sean Brennan
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

%% Demonstration case 1: 2D polytope with 5 edges and 5 vertices
fig_num = 0001;
figure(fig_num);
close(fig_num);

% Fill some vertices with edges. Each vertex is created by its own edge and
% the one prior to it
clear cell_array_edges_in_vertices
Nvertices = 7;
cell_array_edges_in_vertices = cell(Nvertices,1);
for ith_vertex = 1:Nvertices
    raw_edges = [ith_vertex-1; ith_vertex];
    cell_array_edges_in_vertices{ith_vertex,1} = mod(raw_edges,Nvertices)+1;
end

NE = Nvertices;

% Call the function
edge_permutations = fcn_VSkel_findEdgePermutations(cell_array_edges_in_vertices, NE, (fig_num));

% Check variable types
assert(isnumeric(edge_permutations));

% Check lengths
assert(length(edge_permutations(:,1))   <= Nvertices*(Nvertices-2));
assert(length(edge_permutations(1,:))   == 3);


% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));


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


%% Demonstration case: 2D polytope with 2 edges and 2 vertices (produces empty matrix)
fig_num = 1002;
figure(fig_num);
close(fig_num);

% Fill some vertices with edges. Each vertex is created by its own edge and
% the one prior to it
clear cell_array_edges_in_vertices
Nvertices = 2;
cell_array_edges_in_vertices = cell(Nvertices,1);
for ith_vertex = 1:Nvertices
    raw_edges = [ith_vertex-1; ith_vertex];
    cell_array_edges_in_vertices{ith_vertex,1} = mod(raw_edges,Nvertices)+1;
end

NE = Nvertices;

% Call the function
edge_permutations = fcn_VSkel_findEdgePermutations(cell_array_edges_in_vertices, NE, (fig_num));

% Check variable types
assert(isnumeric(edge_permutations));

% Check lengths
assert(isempty(edge_permutations))


% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));

%% Demonstration case: 2D polytope with 3 edges and 3 vertices (gives only 1 permutation)
fig_num = 1003;
figure(fig_num);
close(fig_num);

% Fill some vertices with edges. Each vertex is created by its own edge and
% the one prior to it
clear cell_array_edges_in_vertices
Nvertices = 3;
cell_array_edges_in_vertices = cell(Nvertices,1);
for ith_vertex = 1:Nvertices
    raw_edges = [ith_vertex-1; ith_vertex];
    cell_array_edges_in_vertices{ith_vertex,1} = mod(raw_edges,Nvertices)+1;
end

NE = Nvertices;

% Call the function
edge_permutations = fcn_VSkel_findEdgePermutations(cell_array_edges_in_vertices, NE, (fig_num));

% Check variable types
assert(isnumeric(edge_permutations));

% Check lengths
assert(length(edge_permutations(:,1))   <= Nvertices*(Nvertices-2));
assert(length(edge_permutations(1,:))   == 3);


% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));

%% Demonstration case: 2D polytope with 4 edges and 4 vertices
fig_num = 1004;
figure(fig_num);
close(fig_num);

% Fill some vertices with edges. Each vertex is created by its own edge and
% the one prior to it
clear cell_array_edges_in_vertices
Nvertices = 4;
cell_array_edges_in_vertices = cell(Nvertices,1);
for ith_vertex = 1:Nvertices
    raw_edges = [ith_vertex-1; ith_vertex];
    cell_array_edges_in_vertices{ith_vertex,1} = mod(raw_edges,Nvertices)+1;
end

NE = Nvertices;

% Call the function
edge_permutations = fcn_VSkel_findEdgePermutations(cell_array_edges_in_vertices, NE, (fig_num));

% Check variable types
assert(isnumeric(edge_permutations));

% Check lengths
assert(length(edge_permutations(:,1))   <= Nvertices*(Nvertices-2));
assert(length(edge_permutations(1,:))   == 3);


% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));

%% Demonstration case: 2D polytope with 30 edges and 30 vertices
fig_num = 1005;
figure(fig_num);
close(fig_num);

% Fill some vertices with edges. Each vertex is created by its own edge and
% the one prior to it
clear cell_array_edges_in_vertices
Nvertices = 30;
cell_array_edges_in_vertices = cell(Nvertices,1);
for ith_vertex = 1:Nvertices
    raw_edges = [ith_vertex-1; ith_vertex];
    cell_array_edges_in_vertices{ith_vertex,1} = mod(raw_edges,Nvertices)+1;
end

NE = Nvertices;

% Call the function
edge_permutations = fcn_VSkel_findEdgePermutations(cell_array_edges_in_vertices, NE, (fig_num));

% Check variable types
assert(isnumeric(edge_permutations));

% Check lengths
assert(length(edge_permutations(:,1))   <= Nvertices*(Nvertices-2));
assert(length(edge_permutations(1,:))   == 3);


% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));






%% BREAKING case: 2D polytope with 30 edges and 30 vertices
fig_num = 7001;
figure(fig_num);
close(fig_num);

% Fill some vertices with edges. Each vertex is created by its own edge and
% the one prior to it
clear cell_array_edges_in_vertices

vertices = [0 0; 3 0; 5 5; 7 0; 10 0; 5 10; 0 5; 0 0];

%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTIONALIZE STARTING HERE

if iscell(vertices)
    % Npolytopes = length(vertices);
    %flag_useCells = 1;
else
    % Need to typecast all variables as cell arrays, so methods for cells
    % can be used concurrently whether or not cell array data was given as
    % input.

    Npolytopes = 1;
    vertices = {vertices};
    unit_normal_vectors = {unit_normal_vectors};
    unit_vertex_projection_vectors = {unit_vertex_projection_vectors};
    vector_direction_of_unit_cut = {vector_direction_of_unit_cut};
    flag_vertexIsNonConvex = {flag_vertexIsNonConvex}; 
    max_edge_cuts = {max_edge_cuts};
    flag_useCells = 0;
end

% Is this 2D or 3D?
dimension_of_points = length(vertices{1}(1,:));


% Make a list of all vertices, all vertices in edges, and all edges in
% vertices
all_vertex_positions = []; % The XY(Z) positions of all vertices
all_vertex_polyIDs   = []; % Which polytope each vertex came from
all_vertex_vertexIDs = []; % Which vertex, in the polytope, this vertex came from
cell_array_vertices_in_edges = cell(1,1); % Which vertices are in each edge
cell_array_edges_in_vertices = cell(1,1); % Which edges/faces define a vertex
all_edge_normals = [];

% Create a counting variable to keep track of how many rows were filled by
% previous polytopes, so that rows in current polytope are offset correctly
previous_vertex_offset = 0;

for ith_polytope = 1:Npolytopes

    %%%%%
    % Fill in vertices
    verticesInThisPolytope = vertices{ith_polytope};
    if 2==dimension_of_points
        % If in 2D, need to remove the last point because it's a repeat of
        % the first
        uniqueVerticesThisPolytope = verticesInThisPolytope(1:end-1,:);
    else
        uniqueVerticesThisPolytope = verticesInThisPolytope;
    end    
    all_vertex_positions = [all_vertex_positions; uniqueVerticesThisPolytope]; %#ok<AGROW>

    NuniqueVerticesThisPolytope = length(uniqueVerticesThisPolytope);

    all_vertex_polyIDs = [all_vertex_polyIDs; ones(NuniqueVerticesThisPolytope,1)*ith_polytope]; %#ok<AGROW>
    thisPolyVertexNumbering = (1:NuniqueVerticesThisPolytope)';
    all_vertex_vertexIDs = [all_vertex_vertexIDs; thisPolyVertexNumbering]; %#ok<AGROW>

    %%%%%
    % Fill in edge definitions and edge normals
    if 2==dimension_of_points
        nextEdge = (1:NuniqueVerticesThisPolytope)';
        previousEdge = mod(nextEdge-2,NuniqueVerticesThisPolytope)+1;
        thisVertex = nextEdge;
        nextVertex = mod(thisVertex,NuniqueVerticesThisPolytope)+1;
        for ith_vertex = 1:NuniqueVerticesThisPolytope
            cell_array_edges_in_vertices{ith_vertex + previous_vertex_offset,1} = [previousEdge(ith_vertex,1) nextEdge(ith_vertex,1)]+previous_vertex_offset;
            cell_array_vertices_in_edges{ith_vertex + previous_vertex_offset,1} = [thisVertex(ith_vertex,1) nextVertex(ith_vertex,1)]+previous_vertex_offset;
        end
        all_edge_normals = [all_edge_normals; unit_normal_vectors{ith_polytope}]; %#ok<AGROW>
    else
        error('3D case not yet coded for filling in edges');
    end

    previous_vertex_offset = previous_vertex_offset + NuniqueVerticesThisPolytope;
end % Ends loop through polytopes
%%%%%%%%%%%%%%%%%%%%%

edge_permutations = fcn_VSkel_findEdgePermutations(cell_array_edges_in_vertices, length(cell_array_vertices_in_edges), (fig_num));

Nvertices = length(cell_array_edges_in_vertices);

% Check variable types
assert(isnumeric(edge_permutations));

% Check lengths
assert(length(edge_permutations(:,1))   <= Nvertices*(Nvertices-2));
assert(length(edge_permutations(1,:))   == 3);


% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));



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


%% Demonstration case: 3D polytope (square diamond)
fig_num = 2001;
figure(fig_num);
close(fig_num);

% The following is motivated by a diamond shape in 3D, e.g. a square
% pyramid stacked atop another square pyramid. It has 6 vertices, 4 in the
% middle in a square pattern, and one on each side of the square. It has 8
% external faces, and each vertex connects 4 faces.
clear cell_array_edges_in_vertices
Nvertices = 6;
cell_array_edges_in_vertices{1,1} = [1 2 3 4];
cell_array_edges_in_vertices{2,1} = [1 2 5 6];
cell_array_edges_in_vertices{3,1} = [2 3 6 7];
cell_array_edges_in_vertices{4,1} = [3 4 7 8];
cell_array_edges_in_vertices{5,1} = [4 1 8 5];
cell_array_edges_in_vertices{6,1} = [5 6 7 8];

Nfaces = 8;

% Call the function
edge_permutations = fcn_VSkel_findEdgePermutations(cell_array_edges_in_vertices, Nfaces, (fig_num));

% Check variable types
assert(isnumeric(edge_permutations));

% Check lengths
assert(length(edge_permutations(:,1))   <= (Nfaces-3)*(nchoosek(4,3))*Nvertices);
assert(length(edge_permutations(1,:))   == 4);

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));


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
     
fig_num = 2002;
figure(fig_num);
close(fig_num);

% The following is motivated by a diamond shape in 3D, e.g. a square
% pyramid stacked atop another square pyramid. It has 6 vertices, 4 in the
% middle in a square pattern, and one on each side of the square. It has 8
% external faces, and each vertex connects 4 faces.
clear cell_array_edges_in_vertices
Nvertices = 5;
cell_array_edges_in_vertices{1,1} = [1 2 3 4];
cell_array_edges_in_vertices{2,1} = [1 2 5];
cell_array_edges_in_vertices{3,1} = [2 3 5];
cell_array_edges_in_vertices{4,1} = [3 4 5];
cell_array_edges_in_vertices{5,1} = [4 1 5];

Nfaces = 5;

% Call the function
edge_permutations = fcn_VSkel_findEdgePermutations(cell_array_edges_in_vertices, Nfaces, (fig_num));

% Check variable types
assert(isnumeric(edge_permutations));

% Check lengths
assert(length(edge_permutations(:,1))   <= (Nfaces-3)*(nchoosek(4,3))*Nvertices);
assert(length(edge_permutations(1,:))   == 4);

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));

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
clear cell_array_edges_in_vertices
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
edge_permutations = fcn_VSkel_findEdgePermutations(cell_array_edges_in_vertices, Nfaces, (fig_num));

% Check variable types
assert(isnumeric(edge_permutations));

% Check lengths
assert(length(edge_permutations(:,1))   <= (Nfaces-3)*(nchoosek(4,3))*Nvertices);
assert(length(edge_permutations(1,:))   == 4);



% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));

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
clear cell_array_edges_in_vertices
Nvertices = 4;
cell_array_edges_in_vertices{1,1} = [1 2 3];
cell_array_edges_in_vertices{2,1} = [1 3 4];
cell_array_edges_in_vertices{3,1} = [1 2 4];
cell_array_edges_in_vertices{4,1} = [2 3 4];

Nfaces = 4;

% Call the function
edge_permutations = fcn_VSkel_findEdgePermutations(cell_array_edges_in_vertices, Nfaces, (fig_num));

% Check variable types
assert(isnumeric(edge_permutations));

% Check lengths
assert(length(edge_permutations(:,1))   <= (Nfaces-3)*(nchoosek(4,3))*Nvertices);
assert(length(edge_permutations(1,:))   == 4);



% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));


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
clear cell_array_edges_in_vertices
Nvertices = 30;
cell_array_edges_in_vertices = cell(Nvertices,1);
for ith_vertex = 1:Nvertices
    raw_edges = [ith_vertex-1; ith_vertex];
    cell_array_edges_in_vertices{ith_vertex,1} = mod(raw_edges,Nvertices)+1;
end

NE = Nvertices;

% Call the function
edge_permutations = fcn_VSkel_findEdgePermutations(cell_array_edges_in_vertices, NE, ([]));

% Check variable types
assert(isnumeric(edge_permutations));

% Check lengths
assert(length(edge_permutations(:,1))   <= Nvertices*(Nvertices-2));
assert(length(edge_permutations(1,:))   == 3);

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));

%% Basic example of vertex calculation - non-normal wall shrinking, NO FIGURE, FAST MODE
fig_num = 9902;
figure(fig_num);
close(fig_num);


% Fill some vertices with edges. Each vertex is created by its own edge and
% the one prior to it
clear cell_array_edges_in_vertices
Nvertices = 30;
cell_array_edges_in_vertices = cell(Nvertices,1);
for ith_vertex = 1:Nvertices
    raw_edges = [ith_vertex-1; ith_vertex];
    cell_array_edges_in_vertices{ith_vertex,1} = mod(raw_edges,Nvertices)+1;
end

NE = Nvertices;

% Call the function
edge_permutations = fcn_VSkel_findEdgePermutations(cell_array_edges_in_vertices, NE, (-1));

% Check variable types
assert(isnumeric(edge_permutations));

% Check lengths
assert(length(edge_permutations(:,1))   <= Nvertices*(Nvertices-2));
assert(length(edge_permutations(1,:))   == 3);


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
clear cell_array_edges_in_vertices
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
    edge_permutations = fcn_VSkel_findEdgePermutations(cell_array_edges_in_vertices, NE, ([]));
end
slow_method = toc;

% Do calculation with pre-calculation, FAST_MODE on
tic;
for ith_test = 1:Niterations
    % Call the function
    edge_permutations = fcn_VSkel_findEdgePermutations(cell_array_edges_in_vertices, NE, (-1));
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
