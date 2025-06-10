function [unit_normal_vectors, vector_direction_of_unit_cut, flag_vertexIsNonConvex]  = ...
    fcn_VSkel_polytopeFindUnitDirectionVectors(polytopeStructure,varargin)
%% fcn_VSkel_polytopeFindUnitDirectionVectors
% finds the vector_direction_of_unit_cut to use out of each vertex point,
% e.g. the direction and distance needed to move each point, given a
% unit edge cut
%
% FORMAT:
%
% [unit_normal_vectors, ...
%     vector_direction_of_unit_cut, flag_vertexIsNonConvex] = ...
%     fcn_VSkel_polytopeFindUnitDirectionVectors(polytopeStructure, (fig_num))
%
% INPUTS:
%
%     polytopeStructure: a structure containing the vertices and faces in
%     the MATLAB "patch" format, saved in a subfield called polyPatch, e.g.
%          allVertices = polytopeStructure.polyPatch.Vertices;
%          allFaces    = polytopeStructure.polyPatch.Faces;
%     There are 2 required substructures that contain the above format:
%
%          polytopeStructure.polyPatch: contains a patch representation for
%          all polytopes together, which includes in 3D all faces
%
%          polytopeStructure.subPolyPatch(ith_face): for each face,
%          contains a patch representation using vertices for just that
%          face. In 2D, this corresponds to each object.
%
%    (OPTIONAL INPUTS)
%
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed. As well, if given, this forces the
%      variable types to be displayed as output and as well makes the input
%      check process verbose.
%
% OUTPUTS:
%
%     unit_normal_vectors: an (M+1)-by-2 matrix of the unit vectors that
%     point inward as measured from one vertex to the next. The vector is
%     assumed to be attached to the start of the edge given by vertex M.
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
% DEPENDENCIES:
%
%     fcn_DebugTools_checkInputsToFunctions
%
% % EXAMPLES:
% For additional examples, see: script_test_fcn_VSkel_polytopeFindUnitDirectionVectors
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
% 2025_05_14 by Sean Brennan
% -- added case where calcualtions can include line segments
% 2025_05_14 by Sean Brennan
% -- added case where calcualtions can include one point


% TO DO
% -- none

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
        if nargin < 1 || nargin > 2
            warning('on','backtrace');
            warning('Incorrect number of inputs given: %.0f ',nargin);
            error('Incorrect number of input arguments')
        end

        % Check the polytopeStructure input
        if ~iscell(polytopeStructure) && ~isstruct(polytopeStructure)
            fcn_DebugTools_checkInputsToFunctions(...
                polytopeStructure, '2or3column_of_numbers');
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

if 1==0
    debug_fig_temp = 8575;
    figure(debug_fig_temp);
    clf;

    plot_formatting.vertices_plot.vertexLabels_flagOn = 1;
    plot_formatting.vertices_plot.faceLabels_flagOn   = 1;

    fcn_VSkel_plotPolytopeDetails(...
        polytopeStructure,...
        ([]), ...  % unit_normal_vectors
        ([]), ...  % unit_vertex_projection_vectors
        (plot_formatting),... % plot_formatting
        (debug_fig_temp));  % fig_num

    view(3)
end

% Calculate the faces for each vertex
[allVertices, allEdges, allFaces, facesForEachVertex, edgesTouchingEachVertex_intoVertex, edgesTouchingEachVertex_outofVertex, facesEachEdge] = fcn_INTERNAL_findFacesForEachVertex(polytopeStructure);

% Calculate the unit vectors for each edge/face
[unit_normal_vectors_allFaces, base_points] = fcn_INTERNAL_calcNormalVectors(allVertices, allEdges, allFaces); %#ok<ASGLU>

% Calculate the vector cut directions from unit vector
vector_direction_of_unit_cut_allVertices = fcn_INTERNAL_calcUnitCuts(allVertices, allEdges, allFaces, unit_normal_vectors_allFaces, facesForEachVertex, facesEachEdge);

% Calculate if the vertex is convex
flag_vertexIsNonConvex_allPolys = fcn_INTERNAL_calcNonConvexVertex(edgesTouchingEachVertex_intoVertex, edgesTouchingEachVertex_outofVertex, unit_normal_vectors_allFaces);


% Save results
unit_normal_vectors            = unit_normal_vectors_allFaces;
vector_direction_of_unit_cut   = vector_direction_of_unit_cut_allVertices;
flag_vertexIsNonConvex         = flag_vertexIsNonConvex_allPolys;


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

    plot_formatting.vertices_plot.vertexLabels_flagOn = 1;
    plot_formatting.vertices_plot.faceLabels_flagOn   = 1;
    plot_formatting.vertexIsNonConvex_flagOn          = 1;

    fcn_VSkel_plotPolytopeDetails(...
       polytopeStructure,...
       (unit_normal_vectors), ...  % unit_normal_vectors
       (vector_direction_of_unit_cut), ...  % unit_vertex_projection_vectors
       (plot_formatting),... % plot_formatting
       (fig_num));  % fig_num

    allVertices = polytopeStructure.polyPatch.Vertices;
    dimension_of_points = length(allVertices(1,:));
    if 3==dimension_of_points
        view(3)
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

%% fcn_INTERNAL_findFacesForEachVertex
function [allVertices, allEdges, allFaces, ...
    facesTouchingEachVertex, edgesTouchingEachVertex_intoVertex, edgesTouchingEachVertex_outofVertex, facesEachEdge] = fcn_INTERNAL_findFacesForEachVertex(polytopeStructure)
% Finds which faces each vertex belongs to so that, in subsequent steps,
% one can find the projection vectors. 
% 
% In 2D, the "faces" are the edges that are defined from each vertex.
% In 3D, the faces are the true faces of the 3D object

% Define commonly used variables
allVertices = polytopeStructure.polyPatch.Vertices;
allFaces    = polytopeStructure.polyPatch.Faces;
dimension_of_points = length(allVertices(1,:));
Nvertices = length(allVertices(:,1));

% Fill in all the edges. Each edge is a list of indices as 2 columns, where
% the index in column 1 is the point index the edge starts, and the index
% in column 2 is the point index where the edge ends. 
%
% The faces_each_edge is a listing of the face in which the edge came from.

% In 2D, the "faces" structure is really the edges. We need to build this
% from the saved faces
allEdges = []; % Initialize output
sourceFaceForEdge = []; % Initialize output
for ith_subFace = 1:length(polytopeStructure.subPolyPatch)
    thisFace = polytopeStructure.subPolyPatch(ith_subFace).Faces;
    allEdges = [allEdges; thisFace];        %#ok<AGROW>
    sourceFaceForEdge = [sourceFaceForEdge; ith_subFace*ones(length(thisFace(:,1)),1)]; %#ok<AGROW>
end

facesEachEdge = nan(length(allEdges(:,1)),2);
if 3==dimension_of_points
    % Find matching faces, since each edge has 2
    for ith_edge = 1:length(allEdges(:,1))
        thisEdge = allEdges(ith_edge,:);
        facesEachEdge(ith_edge,1) = sourceFaceForEdge(ith_edge,1);
        edgeBackwards = fliplr(thisEdge);
        [~,rowIndex] = intersect(allEdges,edgeBackwards,'rows');
        if isempty(rowIndex)
            error('Edge found in 3D that has no matching reverse traversal.')
        end
        if length(rowIndex)>2
            error('More than two faces in 3D found that share the same edge.');
        end
        facesEachEdge(ith_edge,2) = sourceFaceForEdge(rowIndex,1);
    
    end
else
    
end

if 2==dimension_of_points
    LongestFace = 2;
    allFaces = allEdges;
    % Initialize outputs
    facesTouchingEachVertex      = nan(Nvertices, LongestFace);
    edgesTouchingEachVertex_intoVertex = nan(Nvertices, 1);
    edgesTouchingEachVertex_outofVertex   = nan(Nvertices, 1);
else
    LongestFace = length(allFaces(1,:));
    % Initialize outputs
    facesTouchingEachVertex      = nan(Nvertices, LongestFace);
    edgesTouchingEachVertex_intoVertex = nan(Nvertices, LongestFace);
    edgesTouchingEachVertex_outofVertex   = nan(Nvertices, LongestFace);
end


% Loop through vertices, finding which faces are involved
for ith_vertex = 1:Nvertices
    % Find all faces that use this vertex
    flag_facesWithThisVertex      = allFaces==ith_vertex;
    flag_edgesWithThisVertex_atStart = allEdges(:,1)==ith_vertex;
    flag_edgesWithThisVertex_atEnd   = allEdges(:,2)==ith_vertex;
    
    summedFlags = sum(flag_facesWithThisVertex,2);
    facesWithThisVertex         = (find(summedFlags))';
    edgesWithThisVertex_atStart = (find(flag_edgesWithThisVertex_atStart))';
    edgesWithThisVertex_atEnd   = (find(flag_edgesWithThisVertex_atEnd))';
    
    if 2==dimension_of_points
        if length(facesWithThisVertex)>2
            warning('more than 2 faces found on a vertex. This should not happen');
        end
        
        facesTouchingEachVertex(ith_vertex,:)      = facesWithThisVertex;
        edgesTouchingEachVertex_intoVertex(ith_vertex,:) = edgesWithThisVertex_atEnd;
        edgesTouchingEachVertex_outofVertex(ith_vertex,:)   = edgesWithThisVertex_atStart;
    else
        NfacesThisVertex = length(facesWithThisVertex);
        facesTouchingEachVertex(ith_vertex,1:NfacesThisVertex)      = facesWithThisVertex;
        NedgesThisVertex = length(edgesWithThisVertex_atEnd);
        if NedgesThisVertex~=length(edgesWithThisVertex_atStart)
            error('Incompatible number of edges into and out of a vertex');
        end
        edgesTouchingEachVertex_intoVertex(ith_vertex,1:NedgesThisVertex) = edgesWithThisVertex_atEnd;
        edgesTouchingEachVertex_outofVertex(ith_vertex,1:NedgesThisVertex)   = edgesWithThisVertex_atStart;

    end


end % Ends for loop

end % Ends fcn_INTERNAL_findFacesForEachVertex


%% fcn_INTERNAL_calcNormalVectors
function [unit_normal_vectors_allFaces, base_points] = fcn_INTERNAL_calcNormalVectors(allVertices, allEdges, allFaces)
% Calculates the unit normal vectors.
% In 2D, this is the unit normal vector to each edge
% In 3D, this is the unit normal vector to the plane


debug_fig = [];

% Is this 2D or 3D?
dimension_of_points = length(allVertices(1,:));

% Is this a 2D (edge) or 3D (face) calculation?
if 2==dimension_of_points
    % Calculate the unit vectors for each edge
    % Get all the indicies for all faces
    startPointIndices = allEdges(:,1);
    endPointIndicies  = allEdges(:,2);
    startVerticies = allVertices(startPointIndices,:);
    endVerticies   = allVertices(endPointIndicies,:);
    base_points = (startVerticies + endVerticies)./2;

    % find distances and unit vectors from vertex to vertex
    difference_vertex_to_vertex = endVerticies - startVerticies;
    distances_vertex_to_vertex = sum(difference_vertex_to_vertex.^2,2).^0.5;
    unit_vectors_vertex_to_vertex = difference_vertex_to_vertex./distances_vertex_to_vertex;

    % Rotate by 90 degrees to get unit normal vectors
    unit_normal_vectors_allFaces = unit_vectors_vertex_to_vertex*[0 1; -1 0];

elseif 3==dimension_of_points
    % Initialize outputs
    Nfaces    = length(allFaces(:,1));
    unit_normal_vectors_allFaces = nan(Nfaces, 3);
    base_points = nan(Nfaces,3);

    % If 3D, have to loop through faces
    % Loop through faces
    for ith_face = 1:Nfaces
        thisFace = allFaces(ith_face,:);
        thisFace = thisFace(~isnan(thisFace));
        thisFaceVertices = allVertices(thisFace,:);

        % Find normal vector for this face
        % [unit_normal_vector, base_point, flags_in_directional_agreement, flags_in_magnitude_agreement] = fcn_geometry_findPlaneNormal(points,(98756));
        [unit_normal_vector, base_point, ~, flags_in_magnitude_agreement] = fcn_geometry_findPlaneNormal(thisFaceVertices,(debug_fig));
        if ~all(flags_in_magnitude_agreement==1)
            error('Face specification given where all the points were not in the same planar face. Unable to proceed.');
        end
        unit_normal_vectors_allFaces(ith_face,:) = unit_normal_vector;
        base_points(ith_face,:) = base_point;



    end
else
    warning('on','backtrace');
    warning('A vector was given that has dimension: %.0d, where 2D was expected',dimension_of_points);
    error('Function not yet coded for anything other than 2D');
end


end % Ends fcn_INTERNAL_calcNormalVectors

%% fcn_INTERNAL_calcUnitCuts
function vector_direction_of_unit_cut_allVertices = fcn_INTERNAL_calcUnitCuts(allVertices, allEdges, allFaces, unit_normal_vectors_allFaces, facesForEachVertex, facesEachEdge)

% Is this 2D or 3D?
dimension_of_points = length(allVertices(1,:));

Nvertices = length(allVertices(:,1));
Nfaces    = length(allFaces(:,1));

% The method of finding unit cuts is different for 2D versus 3D
if 2==dimension_of_points
    % In 2D, the vector projection for a vertex is defined by the 2 edges, one
    % leading into the vertex, and one leading out of the vertex.  

    % Initialize outputs. The vector directions are arranged so that each row
    % is one vertex's direction vector.

    vector_direction_of_unit_cut_allVertices = nan(Nvertices,dimension_of_points);
    for ith_vertex = 1:length(allVertices)

        facesWithThisVertex = (facesForEachVertex(ith_vertex,:))';
        facesWithThisVertex = facesWithThisVertex(~isnan(facesWithThisVertex));

        vector_direction_of_unit_cut_allVertices(ith_vertex,:) = fcn_INTERNAL_calculateUnitCutDirectionGivenFaces(unit_normal_vectors_allFaces,  ith_vertex, facesWithThisVertex, 2);
    end % Ends for loop for 2D vertices

else % This is 3D, so need to loop through it accordingly

    % In 3D, the vector projection for a vertex is defined by 3 planes:
    % 1) the plane defining the current face
    % 2) the plane defining the edge within the current face leading into the vertex
    % 3) the plane defining the edge within the current face leaving the vertex
    % Of note: a vertex can have more than one vertex projection,
    % depending on the face to which it belongs. The maximum number of
    % vertex projections is therefore always less than or equal to the
    % number of faces.
    % Thus, to find the vertex projection, one can loop through faces
    % and, for each vertex in the face, find the edges leading in and
    % out of them. For each edge, find the other faces that share that
    % edge. These define the 3 planes in 3D.

    % Initialize outputs. The vector directions are arranged so that each
    % cell row is one vertex's direction vector.

    vector_direction_of_unit_cut_allVertices = cell(Nvertices,1);

    for ith_face = 1:Nfaces
        thisFace = allFaces(ith_face,:);
        thisFace = thisFace(~isnan(thisFace));
        NverticesThisFace = length(thisFace);
        for ith_vertex = 1:NverticesThisFace

            % Pull out the vertex numbers for the previous, current, and
            % next vertex in this face
            previousIndex = mod(ith_vertex-2,NverticesThisFace)+1;
            previousVertexIndex = thisFace(1,previousIndex);
            currentVertexIndex = thisFace(1,ith_vertex);
            nextIndex = mod(ith_vertex,NverticesThisFace)+1;
            nextVertexIndex = thisFace(1,nextIndex);
            

            % Find all the faces that contain any of these three verticies
            facesWithVertices = 1.0*(allFaces==previousVertexIndex) + 1.0*(allFaces==currentVertexIndex) + 1.0*(allFaces==nextVertexIndex);
            thisVertexFaces = find(sum(facesWithVertices,2)>1);

            % % Find which edges match this sequence
            % previousEdgeRow = [previousVertexIndex currentVertexIndex];
            % nextEdgeRow     = [currentVertexIndex nextVertexIndex];
            % [~,previousEdgeIndex] = intersect(allEdges,previousEdgeRow,'rows');
            % [~,nextEdgeIndex] = intersect(allEdges,nextEdgeRow,'rows');
            % if isempty(previousEdgeIndex) || isempty(nextEdgeIndex)
            %     error('Edge not found');
            % end
            % 
            % % Pull out the faces for previous and next edge
            % facesPrevious = facesEachEdge(previousEdgeIndex,:);
            % facesNext = facesEachEdge(nextEdgeIndex,:);
            % thisVertexFaces= unique([facesPrevious facesNext]);
            if length(thisVertexFaces)<3
                error('Insufficient number of faces found');
            end

            disp('URHERE');

            % Solve for the unit cut direction
            thisDirection = fcn_INTERNAL_calculateUnitCutDirectionGivenFaces(unit_normal_vectors_allFaces, [],  thisVertexFaces, 3);

            currentDirections =  vector_direction_of_unit_cut_allVertices{currentVertexIndex,1};
            currentDirections = [currentDirections; thisDirection] ; %#ok<AGROW>

            vector_direction_of_unit_cut_allVertices{currentVertexIndex,1} = currentDirections;

        end
               
    end % Ends looping through faces in 3D
    
    % Loop through the vectors, keeping only unique ones
    totalNVectors = 0;
    unique_vector_direction_of_unit_cut_allVertices = cell(Nvertices,1);
    for ith_vector = 1:Nvertices
        thisVertexVectors = vector_direction_of_unit_cut_allVertices{ith_vector,1};
        uniqueVectors = unique(thisVertexVectors,'rows','legacy');
        unique_vector_direction_of_unit_cut_allVertices{ith_vector,1} = uniqueVectors;
        totalNVectors = totalNVectors + length(uniqueVectors(:,1));
    end
    vector_direction_of_unit_cut_allVertices = unique_vector_direction_of_unit_cut_allVertices;

    % Convert to matrix?
    if totalNVectors == Nvertices
        temp_vector_direction_of_unit_cut_allVertices = nan(Nvertices,dimension_of_points);
        for ith_vector = 1:Nvertices
            thisVertexVectors = vector_direction_of_unit_cut_allVertices{ith_vector,1};
            temp_vector_direction_of_unit_cut_allVertices(ith_vector,:) = thisVertexVectors;
        end
        vector_direction_of_unit_cut_allVertices = temp_vector_direction_of_unit_cut_allVertices;
    end

end % Ends if statement checking if 2D or 3D



end % Ends fcn_INTERNAL_calcUnitCuts

%% fcn_INTERNAL_calcNonConvexVertex
function flag_vertexIsNonConvex = fcn_INTERNAL_calcNonConvexVertex(edgesTouchingEachVertex_intoVertex, edgesTouchingEachVertex_outofVertex, unit_normal_vectors_allFaces)

% Is this 2D or 3D?
dimension_of_points = length(unit_normal_vectors_allFaces(1,:));
Nvertices = length(edgesTouchingEachVertex_intoVertex(:,1));

if 2==dimension_of_points
    % Find all faces for this vertex
    startFaces = edgesTouchingEachVertex_intoVertex(:,1);
    endFaces   = edgesTouchingEachVertex_outofVertex(:,1);

    % Find their cross product
    cross_products = cross([unit_normal_vectors_allFaces(startFaces,:) zeros(Nvertices,1)],[unit_normal_vectors_allFaces(endFaces,:) zeros(Nvertices,1)],2);
    cross_product_results = cross_products(:,3);

    % Vertex is non convex if cross product less than zero
    flag_vertexIsNonConvex = cross_product_results<0;

else
    % Unclear if nonconvex labels needed for 3D? Solve this later if so.
    flag_vertexIsNonConvex = zeros(Nvertices,1);
end



end % Ends fcn_INTERNAL_calcNonConvexVertex

%% fcn_INTERNAL_calculateUnitCutDirectionGivenFaces
function vector_solution = fcn_INTERNAL_calculateUnitCutDirectionGivenFaces(unit_normal_vectors_allFaces, ith_vertex, facesWithThisVertex, dimension_of_points)
% Set up A matrix for linear equation solution
Amatrix = unit_normal_vectors_allFaces(facesWithThisVertex,:);

if any(isnan(Amatrix),'all') || any(isinf(Amatrix),'all') || rank(Amatrix)<dimension_of_points
    % Check the line segment case. Do the vector directions point
    % in opposite directions? If so, use this direction.
    if isequal(Amatrix(1,:),-1*Amatrix(2,:))
        if 2==dimension_of_points
            vector_solution = unit_normal_vectors_allFaces(ith_vertex,:)*[0 -1; 1 0];
        else
            error('not coded yet');
        end
    else
        % Check the point case - return all nans
        if all(isnan(Amatrix))
            % Do nothing - returns NaN 
            vector_solution = nan(1,dimension_of_points);
        else
            error('Degenerate A matrix found');
        end
    end
else
    % Solve linear equation
    Npoints = length(facesWithThisVertex);
    vector_solution = (Amatrix\ones(Npoints,1))';
end

end % Ends fcn_INTERNAL_calculateUnitCutDirectionGivenFaces