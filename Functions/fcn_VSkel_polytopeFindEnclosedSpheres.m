function [sphereEdgeRadii, definingBoundaries] = ...
    fcn_VSkel_polytopeFindEnclosedSpheres(vertices, unit_normal_vectors, unit_vertex_projection_vectors, vector_direction_of_unit_cut, flag_vertexIsNonConvex, max_edge_cuts, varargin)

%% fcn_VSkel_polytopeFindEnclosedSpheres
% finds the enclosed sphere radii for each vertex point, to each other edge
%
% FORMAT:
%
% [sphereRadii, definingBoundaries] = ...
%     fcn_VSkel_polytopeFindEnclosedSpheres(vertices, unit_normal_vectors, unit_vertex_projection_vectors, flag_vertexIsNonConvex, (fig_num))
%
% INPUTS:
%
%     vertices: a (M+1)-by-2 matrix of xy points with row1 = rowm+1, where
%     M is the number of the individual polytope vertices. The verticies
%     input can also be a cell array of vertex sequences, where each cell
%     array represents a different polytope. If a cell array is given as a
%     vertices input, the outputs are grouped as cell arrays corresponding
%     to the same polytopes. As well, the unit_normal_vectors, 
%     unit_vertex_projection_vectors, etc. are assumed to be cell arrays also.
%
%     unit_normal_vectors: a cell array of dimension M, where
%     each index 1:M stores a N x 2 array of the unit vectors that point
%     inward as measured from one vertex to the next.
%
%     unit_vertex_projection_vectors: a cell array of M, where each index 1:M
%     stores a N x 2 array of the unit vectors that point
%     away from the vertices into the nested shape inside, with M = 1 being
%     the starting unit vectors and N being smaller and smaller for each M value.
%
%     vector_direction_of_unit_cut: an (M+1)-by-2 matrix of the unit
%     vectors that define the magnitude and diretion of the vertices
%     movement into the nested shape inside, assuming a unit magnitude cut.
%     The vector is assumed to be attached to the start of the edge given
%     by vertex M.
%
%     flag_vertexIsNonConvex: an N x 1 array of flags (true or false) that
%     indicate whether the vertex is not convex (1 = NOT convex)
%
%     max_edge_cuts: a (M+1)-by-1 matrix of distances, for each vertex,
%     that can be cut from the edges before the edge "self-intersects",
%     e.g. the projection from each side of the edges meet. If there is not
%     any intersection, the distance is infinite.
%
%     (OPTIONAL INPUTS)
%
%     fig_num: a figure number to plot results. If set to -1, skips any
%     input checking or debugging, no figures will be generated, and sets
%     up code to maximize speed. As well, if given, this forces the
%     variable types to be displayed as output and as well makes the input
%     check process verbose.
%
% OUTPUTS:
%
%     sphereEdgeRadii: a cell array of dimension N containing, in each
%     cell, an array of radii that cause that vertex projection to contact
%     an edge. In each cell array, the are Mx1 radii vectors, where M is =
%     N-2 and N is the number of vertices. The radii are ordered so that
%     the first radii cell array corresponds to the first vertex, etc.
%
%     definingBoundaries: a cell array of dimension N containing, in each cell,
%     an array of which edges constrain each radius of each vertex. In each
%     cell array, there are Mx1 defining edges, where M is = N-2 and N is
%     the number of vertices. The defining edges match the radii ordering,
%     e.g. vertex 2's 3rd sphereRadii edge interaction ID will be in cell
%     array 2, in the 3rd row.
%
% DEPENDENCIES:
%
%     fcn_DebugTools_checkInputsToFunctions
%
% % EXAMPLES:
% For additional examples, see: script_test_fcn_VSkel_polytopeFindEnclosedSpheres
%
% This function was written on 2025_05_02 by S. Brennan
% Questions or comments? sbrennan@psu.edu
%

% Revision History:
% 2022_02_13 - S. Brennan
% -- first write of code
% -- pulled the function out of edge shrinking code
% 2025_04_25 by Sean Brennan
% -- added global debugging options
% -- switched input checking to fcn_DebugTools_checkInputsToFunctions
% 2025_05_02 by Sean Brennan
% -- pulled code out of polytopeFindVertexSkeleton to allow stand-alone
% testings



% TO DO
% -- speed up by checking off which edges were checked previously, and not
% repeating the calculations again for these. For example, in a 4-sided 2D
% polytope, vertex 1 contains edges 41. If this hits edge 3, the result
% will be the same as vertex 4 which is edges 34 hitting edge 1. They all
% contain edges 134.

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==7 && isequal(varargin{end},-1))
    flag_do_debug = 0; % % % % Flag to plot the results for debugging
    flag_check_inputs = 0; % Flag to perform input checking
    flag_max_speed = 1;
else
    % Check to see if we are externally setting debug mode to be "on"
    flag_do_debug = 0; % % % % Flag to plot the results for debugging
    flag_check_inputs = 1; % Flag to perform input checking
    MATLABFLAG_VSKEL_FLAG_CHECK_INPUTS = getenv("MATLABFLAG_VSKEL_FLAG_CHECK_INPUTS");
    MATLABFLAG_VSKEL_FLAG_DO_DEBUG = getenv("MATLABFLAG_VSKEL_FLAG_DO_DEBUG");
    if ~isempty(MATLABFLAG_VSKEL_FLAG_CHECK_INPUTS) && ~isempty(MATLABFLAG_VSKEL_FLAG_DO_DEBUG)
        flag_do_debug = str2double(MATLABFLAG_VSKEL_FLAG_DO_DEBUG);
        flag_check_inputs  = str2double(MATLABFLAG_VSKEL_FLAG_CHECK_INPUTS);
    end
end

% flag_do_debug = 1;

if flag_do_debug
    st = dbstack; %#ok<*UNRCH>
    fprintf(1,'STARTING function: %s, in file: %s\n',st(1).name,st(1).file);
    debug_fig_num = 999978; %#ok<NASGU>
else
    debug_fig_num = []; %#ok<NASGU>
end

%% check input arguments?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____                   _
%  |_   _|                 | |
%    | |  _ __  _ __  _   _| |_ ___
%    | | | '_ \| '_ \| | | | __/ __|
%   _| |_| | | | |_) | |_| | |_\__ \
%  |_____|_| |_| .__/ \__,_|\__|___/
%              | |
%              |_|
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 0==flag_max_speed
    if flag_check_inputs
        % Are there the right number of inputs?
        narginchk(6,7);

        if ~iscell(vertices)
            % Check the vertices input
            fcn_DebugTools_checkInputsToFunctions(...
                vertices, '2or3column_of_numbers');

            NumUniqueVerticies = length(vertices(:,1));

            % Check the unit_normal_vectors input
            fcn_DebugTools_checkInputsToFunctions(...
                unit_normal_vectors, '2or3column_of_numbers',NumUniqueVerticies);

            % Check the unit_vertex_projection_vectors input
            fcn_DebugTools_checkInputsToFunctions(...
                unit_vertex_projection_vectors, '2or3column_of_numbers',NumUniqueVerticies);

            % Check the unit_vertex_projection_vectors input
            fcn_DebugTools_checkInputsToFunctions(...
                vector_direction_of_unit_cut, '2or3column_of_numbers',NumUniqueVerticies);

            % Check the flag_vertexIsNonConvex input
            fcn_DebugTools_checkInputsToFunctions(...
                flag_vertexIsNonConvex*1.00, '1column_of_numbers',NumUniqueVerticies);

            % Check the max_edge_cuts input
            fcn_DebugTools_checkInputsToFunctions(...
                max_edge_cuts, '1column_of_numbers',NumUniqueVerticies);
        end
    end
end

% Does user want to show the plots?
flag_do_plot = 0; % Default is no plotting
if  7 == nargin && (0==flag_max_speed) % Only create a figure if NOT maximizing speed
    temp = varargin{end}; % Last argument is always figure number
    if ~isempty(temp) % Make sure the user is not giving empty input
        fig_num = temp;
        flag_do_plot = 1; % Set flag to do plotting
    end
else
    if flag_do_debug % If in debug mode, do plotting but to an arbitrary figure number
        fig = figure;
        fig_for_debug = fig.Number; %#ok<NASGU>
        flag_do_plot = 1;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _
%  |  \/  |     (_)
%  | \  / | __ _ _ _ __
%  | |\/| |/ _` | | '_ \
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTIONALIZE STARTING HERE

if iscell(vertices)
    Npolytopes = length(vertices);
    flag_useCells = 1;
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
all_vector_direction_of_unit_cut = [];
all_max_edge_cuts = [];

% Create a counting variable to keep track of how many rows were filled by
% previous polytopes, so that rows in current polytope are offset correctly
previous_vertex_offset = 0;

for ith_polytope = 1:Npolytopes

    %%%%%
    % Fill in vertices and vertex projections
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
    all_vector_direction_of_unit_cut = [all_vector_direction_of_unit_cut; vector_direction_of_unit_cut{ith_polytope}]; %#ok<AGROW>

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
        all_max_edge_cuts = [all_max_edge_cuts; max_edge_cuts{ith_polytope}];
    else
        error('3D case not yet coded for filling in edges');
    end

    previous_vertex_offset = previous_vertex_offset + NuniqueVerticesThisPolytope;
end % Ends loop through polytopes
%%%%%%%%%%%%%%%%%%%%%

[edge_permutations, edge_not_in_vertex] = fcn_VSkel_findEdgePermutations(cell_array_edges_in_vertices, length(cell_array_vertices_in_edges), (fig_num));

% For each of the permutations, solve for distances and centers
edge_solutions = zeros(length(edge_permutations(:,1)),3);

for ith_permutation = 1:length(edge_permutations(:,1))
    thisPermutation = edge_permutations(ith_permutation,:);
    vertices_thisPermutation = all_vertex_positions(thisPermutation,:);
    normals_thisPermutation = all_edge_normals(thisPermutation,:);

    % Set up the linear equation, and solve
    y = sum(vertices_thisPermutation.*normals_thisPermutation,2);
    A = [normals_thisPermutation -ones(3,1)];
    if rank(A)>2
        edge_solutions(ith_permutation,:) = (A\y)';
    else
        edge_solutions(ith_permutation,:) = [nan nan nan];
    end    
end

% Remove solutions with negative radii
bad_indicies = find(edge_solutions(:,3)<=0);
if ~isempty(bad_indicies)
    edge_solutions(bad_indicies,:) = nan; 
end

% Remove solutions where radius calculated is larger than the edge will
% project. The minumum cut over the entire column is the minimum allowed.
max_cut_each_permutation = min(all_max_edge_cuts(edge_permutations),[],2); 
bad_indicies = find(edge_solutions(:,3)>(max_cut_each_permutation+eps*1E4));
if ~isempty(bad_indicies)
    edge_solutions(bad_indicies,:) = nan; 
end


% Check if sphereCenter is valid by checking if it is on the edge
% created by projecting each non-vertex edge "inward" by the cut distance
index_of_vertices_of_each_edge = reshape([cell_array_vertices_in_edges{edge_not_in_vertex}], 2,length(edge_not_in_vertex(:,1)))';
edgePointStart = all_vertex_positions(index_of_vertices_of_each_edge(:,1),:);
edgePointEnd   = all_vertex_positions(index_of_vertices_of_each_edge(:,2),:);
edge_vector_direction_of_unit_cut_start = all_vector_direction_of_unit_cut(index_of_vertices_of_each_edge(:,1),:);
edge_vector_direction_of_unit_cut_end   = all_vector_direction_of_unit_cut(index_of_vertices_of_each_edge(:,2),:);

new_edgePointStart = edgePointStart + edge_solutions(:,3).*edge_vector_direction_of_unit_cut_start;
new_edgePointEnd   = edgePointEnd   + edge_solutions(:,3).*edge_vector_direction_of_unit_cut_end;
isOnEdge = fcn_INTERNAL_isPointOnEdge(new_edgePointStart,new_edgePointEnd,edge_solutions(:,1:2), flag_do_debug, debug_fig_num);
bad_indicies = find(isOnEdge==0);
if ~isempty(bad_indicies)
    edge_solutions(bad_indicies,:) = nan;
end


% Initialize outputs
sphereEdgeRadii_allPolytopes    = cell(Npolytopes, 1);
sphereEdgeCenters_allPolytopes  = cell(Npolytopes, 1);
definingBoundaries_allPolytopes = cell(Npolytopes, 1);


for ith_polytope = 1:Npolytopes

    thisPolytopeVertices = vertices{ith_polytope};
    NumUniqueVerticies = length(thisPolytopeVertices(:,1))-1;

    if 2==dimension_of_points

        % Initialize variables
        sphereEdgeRadii_thisPoly      = cell(NumUniqueVerticies+1,1);
        sphereEdgeCenters_thisPoly    = cell(NumUniqueVerticies+1,1);
        definingBoundaries_thisPoly   = cell(NumUniqueVerticies+1,1);

        % For each vertex, solve for the radii to all the non-participating
        % edges. Non-participating edges are those that are not to either side
        % of this vertex.
        for ith_vertex = 1:NumUniqueVerticies

            % Find which edges are intersecting with the current vertex
            [radiiFromVertexToEdge, sphereEdgeCenterArray, edgesConstrainingRadii] = ...
                fcn_INTERNAL_checkVertexProjectionsOntoEdges(thisPolytopeVertices, ith_vertex, unit_normal_vectors{ith_polytope}, unit_vertex_projection_vectors{ith_polytope}, vector_direction_of_unit_cut{ith_polytope}, max_edge_cuts{ith_polytope});

            % sort results so we can compare them next
            [sorted_edgesConstrainingRadii, index_sorted_order] = sort(edgesConstrainingRadii);
            sorted_sphereEdgeCenterArray = sphereEdgeCenterArray(index_sorted_order,:);
            sorted_radiiFromVertexToEdge = radiiFromVertexToEdge(index_sorted_order,:);

            %%%%
            % Compare to linear matrix solution
            % Find edges for this vertex
            edgesThisVertex = cell_array_edges_in_vertices{ith_vertex}; % Which edges are in this vertex?

            % Find matching permutation
            flag_hasFirstIndex = edge_permutations==edgesThisVertex(1,1);
            flag_hasSecondIndex = edge_permutations==edgesThisVertex(1,2);
            flag_permutationsContainEdges = sum(flag_hasFirstIndex,2).*sum(flag_hasSecondIndex,2);
            permutation_indices = find(flag_permutationsContainEdges);
            if isempty(permutation_indices)
                warning('on','backtrace');
                warning('A permutation was not found?. Throwing an error.');
                error('Query permutation not found?');
            end
            constrainingEdge = edge_permutations(permutation_indices,:);
            constrainingEdge(flag_hasFirstIndex(permutation_indices,:)) = 0;
            constrainingEdge(flag_hasSecondIndex(permutation_indices,:)) = 0;
            constrainingEdge = sum(constrainingEdge,2);

            this_solution = edge_solutions(permutation_indices,:);

            fprintf(1,'\n\nVERTEX: %.0d\n',ith_vertex);
            fprintf(1,'Permutation method versus projection method:')            
            table_data = [constrainingEdge, sorted_edgesConstrainingRadii, this_solution(:,3), sorted_radiiFromVertexToEdge, this_solution(:,1), sorted_sphereEdgeCenterArray(:,1), this_solution(:,2), sorted_sphereEdgeCenterArray(:,2)];
            header_strings = [{'Perm Indx'}, {'Proj Indx'},{'Perm R'},{'Proj R'},{'Perm Cx'},{'Proj Cx'},{'Perm Cy'},{'Proj Cy'}]; % Headers for each column
            formatter_strings = [{'%.0d'},{'%.0d'},{'%.3f'},{'%.3f'},{'%.3f'},{'%.3f'},{'%.3f'},{'%.3f'}]; % How should each column be printed?
            N_chars = [10, 10, 10, 10, 10, 10, 10, 10]; % Specify spaces for each column
            fcn_DebugTools_debugPrintTableToNCharacters(table_data, header_strings, formatter_strings,N_chars);

            % Save matricies to cell arrays for this vertex
            sphereEdgeRadii_thisPoly{ith_vertex}     = radiiFromVertexToEdge;
            sphereEdgeCenters_thisPoly{ith_vertex}   = sphereEdgeCenterArray;
            definingBoundaries_thisPoly{ith_vertex}  = edgesConstrainingRadii;

        end

        % Repeat last value as first, since this is the repeated 1st vertex
        
        sphereEdgeRadii_thisPoly{end}     = sphereEdgeRadii_thisPoly{1};
        sphereEdgeCenters_thisPoly{end}   = sphereEdgeCenters_thisPoly{1};
        definingBoundaries_thisPoly{end}  = definingBoundaries_thisPoly{1};

        % Save results into all polytope structure
        sphereEdgeRadii_allPolytopes{ith_polytope}    = sphereEdgeRadii_thisPoly;
        sphereEdgeCenters_allPolytopes{ith_polytope}  = sphereEdgeCenters_thisPoly;
        definingBoundaries_allPolytopes{ith_polytope} = definingBoundaries_thisPoly;

    else
        warning('on','backtrace');
        warning('A vector was given that has dimension: %.0d, where 2D was expected',dimension_of_points);
        error('Function not yet coded for anything other than 2D');
    end
end % Ends loop through polytopes

% Save results
if 0==flag_useCells
    sphereEdgeRadii    = sphereEdgeRadii_allPolytopes{1};
    sphereEdgeCenters  = sphereEdgeCenters_allPolytopes{1};
    definingBoundaries = definingBoundaries_allPolytopes{1};
else
    sphereEdgeRadii    = sphereEdgeRadii_allPolytopes;
    sphereEdgeCenters  = sphereEdgeCenters_allPolytopes;
    definingBoundaries = definingBoundaries_allPolytopes;
end


%% Plot results?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____       _
%  |  __ \     | |
%  | |  | | ___| |__  _   _  __ _
%  | |  | |/ _ \ '_ \| | | |/ _` |
%  | |__| |  __/ |_) | |_| | (_| |
%  |_____/ \___|_.__/ \__,_|\__, |
%                            __/ |
%                           |___/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if flag_do_plot
    
    figure(fig_num);
    clf;
    
    tiledlayout('flow');

    % Find size of vertex domain
    max_XY = max(all_vertex_positions);
    min_XY = min(all_vertex_positions);
    sizePlot = max(max_XY) - min(min_XY);
    nudge = sizePlot*0.006;


    for this_vertex = 1:length(all_vertex_positions(:,1))
        thisPoly   = all_vertex_polyIDs(this_vertex,1);
        thisVertex = all_vertex_vertexIDs(this_vertex,1);

        radiiFromVertexToEdge  = sphereEdgeRadii_allPolytopes{thisPoly}{thisVertex};
        sphereEdgeCenterArray  = sphereEdgeCenters_allPolytopes{thisPoly}{thisVertex};
        edgesConstrainingRadii = definingBoundaries_allPolytopes{thisPoly}{thisVertex};
     
        nexttile;

        % Plot all the polytopes
        for ith_polytope = 1:Npolytopes
            fcn_VSkel_plotPolytopeDetails(...
                vertices{ith_polytope},...
                (unit_normal_vectors{ith_polytope}), ...  % unit_normal_vectors
                ([]), ...  % unit_vertex_projection_vectors
                ([]), ... % vector_direction_of_unit_cut
                (flag_vertexIsNonConvex{ith_polytope}),...  % flag_vertexIsNonConvex
                (1),...  % flag_plotEdgeGhostlines
                (1),...  % flag_plotVertexProjectionGhostlines
                ([]),...  % plot_formatting
                (fig_num));  % fig_num

        end % Ends loop through polytopes


        % Keep the axis from this standard plot, so that all subplots are
        % same
        goodAxis = axis;

        % Plot the current vertex with a red circle
        plot(all_vertex_positions(this_vertex,1), all_vertex_positions(this_vertex,2),'ro','MarkerSize',3);

        % Plot the edge spheres, and label their centers
        for ith_sphere = 1:length(radiiFromVertexToEdge)
            circleCenter = sphereEdgeCenterArray(ith_sphere,:);
            circleRadius = radiiFromVertexToEdge(ith_sphere,1);
            circleEdge   = edgesConstrainingRadii(ith_sphere);

            % Plot the circle center as a large dot, and store the color
            h_fig = plot(circleCenter(1,1),circleCenter(1,2),'.','MarkerSize',20);
            colorUsed = get(h_fig,'Color');

            % Plot the circle boundary in same color
            fcn_geometry_plotCircle(circleCenter,circleRadius,colorUsed,fig_num);

            text(circleCenter(1,1)+nudge, circleCenter(1,2),...
                sprintf('%.0dto %.0d',thisVertex,circleEdge), 'Color',colorUsed);
        end

        % Make axis back to before
        axis(goodAxis);
        title(sprintf('Poly: %.0d, Vertex: %.0d', thisPoly, thisVertex));
    end
end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file);
end

end % Ends INTERNAL_fcn_findUnitDirectionVectors



%% Functions follow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   ______                _   _
%  |  ____|              | | (_)
%  | |__ _   _ _ __   ___| |_ _  ___  _ __  ___
%  |  __| | | | '_ \ / __| __| |/ _ \| '_ \/ __|
%  | |  | |_| | | | | (__| |_| | (_) | | | \__ \
%  |_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|___/
%
% See: https://patorjk.com/software/taag/#p=display&f=Big&t=Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%§

%% fcn_INTERNAL_checkVertexProjectionsOntoEdges
function [radiiFromVertexToEdge, sphereCenterArray, edgesConstrainingRadii] = ...
            fcn_INTERNAL_checkVertexProjectionsOntoEdges(vertices, ith_vertex, unit_normal_vectors, unit_vertex_projection_vectors, vector_direction_of_unit_cut, max_edge_cuts)

% This function loops through edges and checks to see if the test vertex,
% given by vertexPoint, would ever be tangent to each edge. 

flag_do_debug = 0;
debug_fig_num = 67676;

vertexPoint = vertices(ith_vertex,:);
vertexNormal = unit_normal_vectors(ith_vertex,:);
vertexProjection = unit_vertex_projection_vectors(ith_vertex,:);
NumUniqueVerticies = length(vertices(:,1))-1;

% Initialize arrays for outputs
radiiFromVertexToEdge  = zeros(NumUniqueVerticies-2,1);
edgesConstrainingRadii = zeros(NumUniqueVerticies-2,1);
sphereCenterArray      = zeros(NumUniqueVerticies-2,2);

% The edges tested start with the edge ahead of this vertex's edge (which is the one
% to the right of the vertex), and move up until the one before
% this vertex, to the left. This will omit 2 verticies. Including
% the 2 omitted, and the one we are on, we go up by
% NumUniqueVerticies - 3. We convert number back into real edge
% numbering. For example, if there are 4 edges, edge 5 is really
% edge 1, edge 6 is really edge 2, etc. The mod function winds the
% for-loop numbering back into real edge numbering.
edgesToTest = (ith_vertex:(ith_vertex+NumUniqueVerticies-3))';
edgesToTest = mod(edgesToTest,NumUniqueVerticies)+1;

% % What is the maximum vertex distance? This is the most the vertex can move
% % before its edge disappears. It will be the minimum of all the edges to
% % which it belongs.
% previousEdge = mod(ith_vertex-2,NumUniqueVerticies)+1;
% maxAllowableVertexDistance = min(max_edge_cuts(ith_vertex),max_edge_cuts(previousEdge));

% Loop through non-participating edges.
for ith_edge = 1:length(edgesToTest)
    current_edge   = edgesToTest(ith_edge,1);
    edgePointStart = vertices(current_edge,:);
    edgePointEnd   = vertices(current_edge+1,:);
    edgeNormal     = unit_normal_vectors(current_edge,:);
    edge_vector_direction_of_unit_cut_start = vector_direction_of_unit_cut(current_edge,:);
    edge_vector_direction_of_unit_cut_end   = vector_direction_of_unit_cut(current_edge+1,:);

    % Keep the minimum of the current vertex allowable cut distance, and the
    % test edge's minimum distance.
    % maxAllowableDistance = min(maxAllowableVertexDistance,max_edge_cuts(current_edge));
    maxAllowableDistance = max_edge_cuts(current_edge);


    % For debugging
    if 1==flag_do_debug
        figure(debug_fig_num);
        clf;

        fcn_VSkel_plotPolytopeDetails(...
            vertices,...
            (unit_normal_vectors), ...  % unit_normal_vectors
            (unit_vertex_projection_vectors), ...  % unit_vertex_projection_vectors
            ([]), ... % vector_direction_of_unit_cut
            ([]),...  % flag_vertexIsNonConvex
            (1),...  % flag_plotEdgeGhostlines
            (1),...  % flag_plotVertexProjectionGhostlines
            ([]),...  % plot_formatting
            (debug_fig_num));  % fig_num

        goodAxis = axis;

        title(sprintf('Vertex: %.0d, Edge: %.0d',ith_vertex, current_edge));

        % Plot the current vertex with a red circle
        plot(vertices(ith_vertex,1),vertices(ith_vertex,2),'ro','Linewidth',2, 'MarkerSize',10);     

        % Plot the edge being tested in some color, and save the color for
        % later
        edgeData = [edgePointStart; edgePointEnd];
        h_fig = plot(edgeData(:,1),edgeData(:,2),'-','Linewidth',5);
        colorUsed = get(h_fig,'Color');


    end


    % See derivation in PPT documentation for the equations below.
    % Requires diagrams to explain where they come from
    denominator = sum(((vertexNormal - edgeNormal).*vertexProjection),2);
    if abs(denominator)<1E-10
        d = nan;
    else
        d = sum((vertexPoint - edgePointStart).*(edgeNormal),2) / denominator;
    end

    r = d*sum(vertexProjection.*vertexNormal,2);

    % If cut size (radius) is greater than max, set to NaN
    if r>(maxAllowableDistance+eps*1E4) || r<=(0+eps*1E4)
        d = nan;
        r = nan;
    end

    % The sphere's center will be the projection from the vertex
    % point by distance d
    sphereCenter = vertexPoint + vertexProjection*d;

    % Do we need to check if the radius (e.g. the cut distance) is valid?
    if ~isnan(r)
        % Check to see if point is on edge by projecting it back toward
        % the edge by the radius times the edgeNormal
        contactPoint = sphereCenter - edgeNormal*r;


        % For debugging
        if 1==flag_do_debug
            figure(debug_fig_num);

            % Plot the circle center as a large dot, and store the color
            plot(sphereCenter(1,1),sphereCenter(1,2),'.','MarkerSize',20,'Color',colorUsed);

            % Plot the circle boundary in same color
            fcn_geometry_plotCircle(sphereCenter,r,colorUsed,debug_fig_num);

            % Plot the contact point
            plot(contactPoint(1,1),contactPoint(1,2),'.','MarkerSize',30, 'Color',colorUsed);
        end


        % Check if sphereCenter is valid by checking if it is on the edge
        % created by projecting the edge "inward" by the cut distance
        new_edgePointStart = edgePointStart + r*edge_vector_direction_of_unit_cut_start;
        new_edgePointEnd   = edgePointEnd   + r*edge_vector_direction_of_unit_cut_end;

        isOnEdge = fcn_INTERNAL_isPointOnEdge(new_edgePointStart,new_edgePointEnd,sphereCenter, flag_do_debug, debug_fig_num);
        if 1==flag_do_debug
            figure(debug_fig_num);
            axis(goodAxis);
            title(sprintf('Vertex: %.0d, Edge: %.0d, isOnEdge: %.0f',ith_vertex, current_edge,  isOnEdge));
        end
        if ~isOnEdge
            r = nan;
            sphereCenter = [nan nan];
        end

    else
        isOnEdge = 0;
    end

    if 1==flag_do_debug
        figure(debug_fig_num);
        title(sprintf('Vertex: %.0d, Edge: %.0d, isOnEdge: %.0f',ith_vertex, current_edge,  isOnEdge));
    end

    % Save results to matricies
    radiiFromVertexToEdge(ith_edge,1)  = r;
    sphereCenterArray(ith_edge,:) = sphereCenter;
    edgesConstrainingRadii(ith_edge,1) = current_edge;
end

end % Ends fcn_INTERNAL_checkVertexProjectionsOntoEdges

%% fcn_INTERNAL_isPointOnEdge
function isOnEdge = fcn_INTERNAL_isPointOnEdge(edgeStart,edgeEnd,testPoint, flag_do_debug, debug_fig_num)
NtestPoints = length(testPoint(:,1));
isOnEdge = ones(NtestPoints,1);

methodToCheck = 1;
tolerance = 1E-5;

edgeLength     = real(sum((edgeEnd - edgeStart).^2,2).^0.5);

% For debugging

if 1==flag_do_debug && 1==NtestPoints
    figure(debug_fig_num);

    inputdata = [edgeStart; edgeEnd];
    h_fig = plot(inputdata(:,1),inputdata(:,2),'.-','MarkerSize',10,'LineWidth',3);
    colorUsed = get(h_fig,'Color');
end

if methodToCheck==1
    % Use vector projection
    edgeVector = edgeEnd - edgeStart;
    testVector = testPoint - edgeStart;
    unit_edgeVector = edgeVector./edgeLength;


    % Do dot product to see if within length
    dotProduct = sum(unit_edgeVector.*testVector,2);
    flags_notInLength =  (dotProduct<(0-tolerance)) + (dotProduct>(edgeLength+tolerance));
    isOnEdge(find(flags_notInLength>0)) = 0; %#ok<FNDSB>

    % Do dot product to see if not orthogonal
    ortho_edgeVector = unit_edgeVector*[0 1; -1 0];
    dotProduct = sum(ortho_edgeVector.*testVector,2);
    flags_notOrthogonal = (dotProduct<(0-tolerance)) + (dotProduct>(0+tolerance));
    isOnEdge(find(flags_notOrthogonal>0)) = 0; %#ok<FNDSB>
else

    % Use the triangle inequality to check if test point is on edge
    pointDistances = sum((testPoint - edgeStart).^2,2).^0.5 + sum((testPoint - edgeEnd).^2,2).^0.5;
    offsetError = abs(edgeLength - pointDistances);
    if offsetError>tolerance
        isOnEdge = 0;
    else
        isOnEdge = 1;
    end
end

if 1==flag_do_debug && 1==NtestPoints
    figure(debug_fig_num);
    if 1==isOnEdge
        plot(testPoint(:,1),testPoint(:,2),'o','MarkerSize',20,'Color',colorUsed);
    else
        plot(testPoint(:,1),testPoint(:,2),'x','MarkerSize',20,'Color',[1 0 0]);
    end
end

end % Ends fcn_INTERNAL_isPointOnEdge

