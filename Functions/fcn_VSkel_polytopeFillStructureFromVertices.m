function polytopeStructure = fcn_VSkel_polytopeFillStructureFromVertices(vertices, varargin)

%% fcn_VSkel_polytopeFillStructureFromVertices
% given a set of vertices, fills key details for a 2D or 3D polytope in a
% polytopeStructure format, where the subfields of polytopeStructure are
% compatible with MATLAB "patch" formats.
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
%     polytopeStructure: a structure containing the vertices and faces in
%     the MATLAB "patch" format, saved in a subfield called polyPatch, e.g.
%          allVertices = polytopeStructure.polyPatch.Vertices;
%          allFaces    = polytopeStructure.polyPatch.Faces;
%     There are 2 required substructures that contain the above format:
%
%          polytopeStructure.polyPatch: contains a patch representation for
%          all polytopes together, which includes in 3D all faces. The
%          vertices and faces of polyPatch represent the entire object (in
%          3D ) or object field (in 2D)
%
%          polytopeStructure.subPolyPatch(ith_face): for each face,
%          contains a patch representation using vertices for just that
%          face. In 2D, this corresponds to each object. The vertices and
%          faces correspond, respectively, to the vertices of only one face
%          and the edges that define that face.
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
% 2025_05_02 - S. Brennan
% -- first write of code
% -- pulled code out of polytopeFindVertexSkeleton to allow stand-alone
%    testings
% 2025_06_01 - S. Brennan
% -- updated docstrings


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



% Initialize outputs
polytopeStructure = struct;
polyPatch = struct;


% Initialize a list of all vertices, all vertices in edges, and all edges
% in vertices


cell_array_unique_vertices_in_faces = cell(Npolytopes,1); % Which vertices are in each edge
longestFace = 0; % How many columns are needed to represent a face?


all_vertex_positions = []; % Lists all the vertices given by the user, including repeats
all_vertex_faceIDs   = []; % For each user-given vertex, lists which face it was used in
all_vertex_vertexIDs = []; % For each vertex, lists the vertex numbering in the original face



%%%%%
% Loop through faces, and for each, fill in vertices, finding for each
% added face which vertices are unique
for ith_face = 1:Npolytopes
 
    verticesInThisPolytope = vertices{ith_face};
    NuniqueVerticesThisPolytope = length(verticesInThisPolytope);

    if ith_face ==1
        unique_vertex_positions = verticesInThisPolytope; % The XY(Z) positions of all vertices
    end

    % Find which vertices are new
    [new_vertices,~] = setdiff(verticesInThisPolytope,unique_vertex_positions,'rows','stable');
    unique_vertex_positions = [unique_vertex_positions; new_vertices]; %#ok<AGROW>
    all_vertex_positions = [all_vertex_positions; verticesInThisPolytope]; %#ok<AGROW>

    all_vertex_faceIDs = [all_vertex_faceIDs; ones(NuniqueVerticesThisPolytope,1)*ith_face]; %#ok<AGROW>

    [~, ~, thisPolyVertexNumbering] = intersect(verticesInThisPolytope, unique_vertex_positions,'rows','stable');
    all_vertex_vertexIDs = [all_vertex_vertexIDs; thisPolyVertexNumbering]; %#ok<AGROW>

    % How many vertices are in this face?
    longestFace = max(longestFace,length(thisPolyVertexNumbering));

    % Fill in faces
    cell_array_unique_vertices_in_faces{ith_face} = thisPolyVertexNumbering;
end


% Transfer cell array to matrix form, which is required for the Faces field
array_vertices_in_faces = nan(Npolytopes,longestFace); % Which vertices define a face. Each row is a different face. Each row contains the vertices that define the face, in order
for ith_cell = 1:Npolytopes
    thisFaceVertices = cell_array_unique_vertices_in_faces{ith_cell}';
    Nvertices = length(thisFaceVertices);
    array_vertices_in_faces(ith_cell, 1:Nvertices) = thisFaceVertices; % Which vertices define a face
end

%%%%%
% Check if any vertices from one face are embedded within another face
for ith_vertex = 1:length(unique_vertex_positions)
    indicesFaceContainsVertex = array_vertices_in_faces == ith_vertex;
    flag_facesWithoutThisVertex = sum(indicesFaceContainsVertex,2)==0;
    facesToTest = find(flag_facesWithoutThisVertex);
    for ith_test = 1:length(facesToTest)
        faceToTest = facesToTest(ith_test);
        URHERE
        flag_containsVertex = fcn_INTERNAL_checkFaceContainsVertex(array_vertices_in_faces, unique_vertex_positions, ith_vertex, faceToTest);
    end
end



% %%%
% % Find which faces touch each vertex
% NtotalVertices = length(unique_vertex_positions(:,1));
% array_faces_in_vertices = nan(NtotalVertices,Npolytopes); % Which faces touch a vertex. Each row is a vertex. Each column is one of the faces. 
% for ith_vertex = 1:NtotalVertices
%     % Find the vertex in the list
%     faces_flagged = sum(array_vertices_in_faces==ith_vertex,2)>0;
%     array_faces_in_vertices(ith_vertex,faces_flagged) = 1;
% 
% end


% Get the color ordering
colorOrdering = lines; %colormap('parula');

%%%%%%%
% Fill in the structure information for polytopeStructure.polyPatch
polyPatch.Vertices = unique_vertex_positions;
polyPatch.Faces = array_vertices_in_faces;

% Set face coloring
polyPatch.FaceColor = 'flat';
all_faceNumbers = (1:Npolytopes)';
colorIndicesAllFaces = mod(all_faceNumbers-1,length(colorOrdering(:,1)))+1;
polyPatch.FaceVertexCData = colorOrdering(colorIndicesAllFaces,:);
polyPatch.FaceAlpha = 0.3;

% Set edge line style
% polyPatch.LineWidth = 2;

polytopeStructure.polyPatch = polyPatch;

%%%%
% Fill in patch information for each
% polytopeStructure.subPolyPatch(ith_face). This fills in edge definitions
% for each face

for ith_face = 1:Npolytopes

    indicesUsedThisFace_withNaNs = array_vertices_in_faces(ith_face,:);
    indicesUsedThisFace = indicesUsedThisFace_withNaNs(~isnan(indicesUsedThisFace_withNaNs));

    NuniqueVerticesThisPolytope = length(indicesUsedThisFace);
    ordinatesThisFace = (1:NuniqueVerticesThisPolytope)';
    nextOrdinatesThisFace = [ordinatesThisFace(2:end); ordinatesThisFace(1)];

    % Convert ordinates to indices
    vertexIndicesThisFace = indicesUsedThisFace(ordinatesThisFace);
    vertexIndicesNextThisFace = indicesUsedThisFace(nextOrdinatesThisFace);

    facesThisVertex = [vertexIndicesThisFace' vertexIndicesNextThisFace'];
    
    % Save faces for all edges in this polytope
    polyPatch.Vertices = unique_vertex_positions;
    polyPatch.Faces = facesThisVertex;    

    % Set face colorings
    thisFaceColorIndex = colorIndicesAllFaces(ith_face);
    polyPatch.FaceColor = colorOrdering(thisFaceColorIndex,:);
    polyPatch.FaceAlpha = 0.5;
    polyPatch.FaceVertexCData = colorOrdering(ordinatesThisFace,:);

    % Set edge properties
    %polyPatch.EdgeColor = 'flat';
    polyPatch.LineWidth = 1;
    
    polytopeStructure.subPolyPatch(ith_face) = polyPatch;

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

    nexttile;

    % Plot all the polytopes with shading   
    fcn_VSkel_plotPolytopeDetails(...
        polytopeStructure,...
        ([]), ...  % unit_normal_vectors
        ([]), ...  % unit_vertex_projection_vectors
        ([]),...   % plot_formatting
        (fig_num));  % fig_num

    if 3==dimension_of_points
        view(3)
    end

    goodAxis = axis;
    title('All faces');

    for ith_face = 1:Npolytopes

        nexttile;

        clear tempPatch
        tempPatch = struct;
        tempPatch.Vertices = unique_vertex_positions;
        tempPatch.Faces = array_vertices_in_faces(ith_face,:);

        % Set face coloring
        tempPatch.FaceColor = 'flat';
        tempPatch.FaceVertexCData = colorOrdering(ith_face,:);
        tempPatch.FaceAlpha = 0.3;


        patch(tempPatch);
        title(sprintf('Face: %.0d',ith_face));

        xlabel('X');
        ylabel('Y');
        zlabel('Z');


        % Make axis back to before
        axis(goodAxis);

        if 3==dimension_of_points
            view(3)
        end

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
