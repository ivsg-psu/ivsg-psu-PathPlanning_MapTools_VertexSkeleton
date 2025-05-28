function polytopeStructure = fcn_VSkel_polytopeFillStructureFromVertices(vertices, varargin)

%% fcn_VSkel_polytopeFillStructureFromVertices
% given a set of vertices, fills key details for a 2D or 3D polytope. 
%
% FORMAT:
%
% polytopeStructure = fcn_VSkel_polytopeFillStructureFromVertices(vertices,  (fig_num))
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
% For additional examples, see: script_test_fcn_VSkel_polytopeFillStructureFromVertices
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
if (nargin==2 && isequal(varargin{end},-1))
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
        narginchk(1,2);

        if ~iscell(vertices)
            % Check the vertices input
            fcn_DebugTools_checkInputsToFunctions(...
                vertices, '2or3column_of_numbers');

            % NumUniqueVerticies = length(vertices(:,1));
            % 
            % % Check the unit_normal_vectors input
            % fcn_DebugTools_checkInputsToFunctions(...
            %     unit_normal_vectors, '2or3column_of_numbers',NumUniqueVerticies);
            % 
            % % Check the unit_vertex_projection_vectors input
            % fcn_DebugTools_checkInputsToFunctions(...
            %     unit_vertex_projection_vectors, '2or3column_of_numbers',NumUniqueVerticies);
            % 
            % % Check the unit_vertex_projection_vectors input
            % fcn_DebugTools_checkInputsToFunctions(...
            %     vector_direction_of_unit_cut, '2or3column_of_numbers',NumUniqueVerticies);
            % 
            % % Check the flag_vertexIsNonConvex input
            % fcn_DebugTools_checkInputsToFunctions(...
            %     flag_vertexIsNonConvex*1.00, '1column_of_numbers',NumUniqueVerticies);
            % 
            % % Check the max_edge_cuts input
            % fcn_DebugTools_checkInputsToFunctions(...
            %     max_edge_cuts, '1column_of_numbers',NumUniqueVerticies);
        end
    end
end

% Does user want to show the plots?
flag_do_plot = 0; % Default is no plotting
if  2 == nargin && (0==flag_max_speed) % Only create a figure if NOT maximizing speed
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

if iscell(vertices)
    Npolytopes = length(vertices);
else
    % Need to typecast all variables as cell arrays, so methods for cells
    % can be used concurrently whether or not cell array data was given as
    % input.

    Npolytopes = 1;
    vertices = {vertices};
end

% Is this 2D or 3D?
dimension_of_points = length(vertices{1}(1,:));


% Make a list of all vertices, all vertices in edges, and all edges in
% vertices
all_vertex_positions = []; % The XY(Z) positions of all vertices
all_vertex_polyIDs   = []; % Which polytope each vertex came from
all_vertex_vertexIDs = []; % Which vertex, in the polytope, this vertex came from
cell_array_vertices_in_faces = cell(1,1); % Which vertices are in each edge
cell_array_faces_in_vertices = cell(1,1); % Which edges/faces define a vertex
% all_face_normals = [];
% all_vector_direction_of_unit_cut = [];
% all_max_face_cuts = [];

% Initialize outputs
polytopeStructure = struct;
polyPatch = struct;

% Create a counting variable to keep track of how many rows were filled by
% previous polytopes, so that rows in current polytope are offset correctly
previous_vertex_offset = 0;

longestFace = 2; % How many columns are needed to represent a face?

for ith_polytope = 1:Npolytopes
    

    %%%%%
    % Fill in vertices and vertex projections
    verticesInThisPolytope = vertices{ith_polytope};
    all_vertex_positions = [all_vertex_positions; verticesInThisPolytope]; %#ok<AGROW>
    NuniqueVerticesThisPolytope = length(verticesInThisPolytope);

    all_vertex_polyIDs = [all_vertex_polyIDs; ones(NuniqueVerticesThisPolytope,1)*ith_polytope]; %#ok<AGROW>
    thisPolyVertexNumbering = (1:NuniqueVerticesThisPolytope)';
    all_vertex_vertexIDs = [all_vertex_vertexIDs; thisPolyVertexNumbering]; %#ok<AGROW>


    %%%%%
    % Fill in face definitions
    if 2==dimension_of_points
        nextFace = (1:NuniqueVerticesThisPolytope)';
        previousFace = mod(nextFace-2,NuniqueVerticesThisPolytope)+1;
        thisVertex = nextFace;
        nextVertex = mod(thisVertex,NuniqueVerticesThisPolytope)+1;
        for ith_vertex = 1:NuniqueVerticesThisPolytope
            cell_array_faces_in_vertices{ith_vertex + previous_vertex_offset,1} = [previousFace(ith_vertex,1) nextFace(ith_vertex,1)]+previous_vertex_offset;
            cell_array_vertices_in_faces{ith_vertex + previous_vertex_offset,1} = [thisVertex(ith_vertex,1) nextVertex(ith_vertex,1)]+previous_vertex_offset;
        end
    else
        cell_array_vertices_in_faces{ith_polytope,1} = thisPolyVertexNumbering+previous_vertex_offset;        
    end

    %%%%
    % Fill in vectors
    % all_face_normals = [all_face_normals; unit_normal_vectors{ith_polytope}]; %#ok<AGROW>
    % all_max_face_cuts = [all_max_face_cuts; max_edge_cuts{ith_polytope}];

    previous_vertex_offset = previous_vertex_offset + NuniqueVerticesThisPolytope;

end % Ends loop through polytopes

% Get the color ordering
colorOrdering = lines; %colormap('parula');

%%%%%%%
% Fill in the patch details for all edges

% Save faces for all edges
faces = nan(length(cell_array_vertices_in_faces),longestFace);
for ith_face = 1:length(cell_array_faces_in_vertices)
    thisFaceDefinition = cell_array_faces_in_vertices{ith_face};
    faces(ith_face,1:length(thisFaceDefinition)) = thisFaceDefinition;
end

polyPatch.Vertices = all_vertex_positions;
polyPatch.Faces = faces;
% poly.FaceVertexCData = [0; 1; 0.5]; % Fix this later
polyPatch.FaceColor = 'none'; % flat
polyPatch.EdgeColor = 'flat';
% colorIndices = mod(vertexNumbering-1,length(colorOrdering(:,1)))+1;
%vertexNumbering = (1:length(polyPatch.Vertices(:,1)))';
colorIndices = mod(all_vertex_polyIDs-1,length(colorOrdering(:,1)))+1;
polyPatch.FaceVertexCData = colorOrdering(colorIndices,:);
polyPatch.LineWidth = 2;
polytopeStructure.polyPatch = polyPatch;

%%%%
% Fill in patch information for each polytope
for ith_polytope = 1:Npolytopes
    thisPolytopeIndicies = find(all_vertex_polyIDs==ith_polytope);

    % Save faces for all edges in this polytope
    faces = thisPolytopeIndicies';
    polyPatch.Vertices = all_vertex_positions;
    polyPatch.Faces = faces;    
    polyPatch.FaceColor = 'flat';
    polyPatch.FaceAlpha = 0.5;
    polyPatch.EdgeColor = 'flat';
    % colorIndices = mod(vertexNumbering-1,length(colorOrdering(:,1)))+1;
    %vertexNumbering = (1:length(polyPatch.Vertices(:,1)))';
    %colorIndices = mod(all_vertex_polyIDs-1,length(colorOrdering(:,1)))+1;
    polyPatch.FaceVertexCData = colorOrdering(ith_polytope,:);
    polyPatch.LineWidth = 1;
    polytopeStructure.subPolyPatch(ith_polytope) = polyPatch;


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
    fcn_VSkel_plotPolytopeDetails(...
        polytopeStructure,...
        ([]), ...  % unit_normal_vectors
        ([]), ...  % unit_vertex_projection_vectors
        ([]), ... % vector_direction_of_unit_cut
        ([]),...  % flag_vertexIsNonConvex
        (1),...  % flag_plotEdgeGhostlines
        (1),...  % flag_plotVertexProjectionGhostlines
        ([]),...  % plot_formatting
        (fig_num));  % fig_num


    figure(fig_num);
    patch(polytopeStructure.polyPatch);   

    % Plot all the polytopes with shading   
    for ith_polytope = 1:Npolytopes
        patch(polytopeStructure.subPolyPatch(ith_polytope));
    end

    if 3==dimension_of_points
        view(3)
    end

    % Plot all the polytopes    
    % fcn_VSkel_plotPolytopeDetails(poly, (fig_num));  

    % Nothing to do
    % figure(fig_num);
    % clf;
    % 
    % tiledlayout('flow');
    % 
    % % Find size of vertex domain
    % max_XY = max(all_vertex_positions);
    % min_XY = min(all_vertex_positions);
    % sizePlot = max(max_XY) - min(min_XY);
    % nudge = sizePlot*0.006;
    % 
    % 
    % for this_vertex = 1:length(all_vertex_positions(:,1))
    %     thisPoly   = all_vertex_polyIDs(this_vertex,1);
    %     thisVertex = all_vertex_vertexIDs(this_vertex,1);
    % 
    %     radiiFromVertexToEdge  = sphereEdgeRadii_allPolytopes{thisPoly}{thisVertex};
    %     sphereEdgeCenterArray  = sphereEdgeCenters_allPolytopes{thisPoly}{thisVertex};
    %     edgesConstrainingRadii = definingBoundaries_allPolytopes{thisPoly}{thisVertex};
    % 
    %     nexttile;
    % 
    %     % Plot all the polytopes
    %     for ith_polytope = 1:Npolytopes
    %         fcn_VSkel_plotPolytopeDetails(...
    %             vertices{ith_polytope},...
    %             (unit_normal_vectors{ith_polytope}), ...  % unit_normal_vectors
    %             ([]), ...  % unit_vertex_projection_vectors
    %             ([]), ... % vector_direction_of_unit_cut
    %             (flag_vertexIsNonConvex{ith_polytope}),...  % flag_vertexIsNonConvex
    %             (1),...  % flag_plotEdgeGhostlines
    %             (1),...  % flag_plotVertexProjectionGhostlines
    %             ([]),...  % plot_formatting
    %             (fig_num));  % fig_num
    % 
    %     end % Ends loop through polytopes
    % 
    % 
    %     % Keep the axis from this standard plot, so that all subplots are
    %     % same
    %     goodAxis = axis;
    % 
    %     % Plot the current vertex with a red circle
    %     plot(all_vertex_positions(this_vertex,1), all_vertex_positions(this_vertex,2),'ro','MarkerSize',3);
    % 
    %     % Plot the edge spheres, and label their centers
    %     for ith_sphere = 1:length(radiiFromVertexToEdge)
    %         circleCenter = sphereEdgeCenterArray(ith_sphere,:);
    %         circleRadius = radiiFromVertexToEdge(ith_sphere,1);
    %         circleEdge   = edgesConstrainingRadii(ith_sphere);
    % 
    %         % Plot the circle center as a large dot, and store the color
    %         h_fig = plot(circleCenter(1,1),circleCenter(1,2),'.','MarkerSize',20);
    %         colorUsed = get(h_fig,'Color');
    % 
    %         % Plot the circle boundary in same color
    %         fcn_geometry_plotCircle(circleCenter,circleRadius,colorUsed,fig_num);
    % 
    %         text(circleCenter(1,1)+nudge, circleCenter(1,2),...
    %             sprintf('%.0dto %.0d',thisVertex,circleEdge), 'Color',colorUsed);
    %     end
    % 
    %     % Make axis back to before
    %     axis(goodAxis);
    %     title(sprintf('Poly: %.0d, Vertex: %.0d', thisPoly, thisVertex));
    % end
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
