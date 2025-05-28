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


allVertices = polytopeStructure.polyPatch.Vertices;
allFaces    = polytopeStructure.polyPatch.Faces;

% Calculate the faces for each vertex
[facesForEachVertex, facesForEachVertex_from, facesForEachVertex_to] = fcn_INTERNAL_findFacesForEachVertex(allVertices, allFaces);

% Calculate the unit vectors for each edge
unit_normal_vectors_allFaces = fcn_INTERNAL_calcNormalVectors(allVertices, allFaces);

% Calculate the vector directions from unit vector
vector_direction_of_unit_cut_allVertices = fcn_INTERNAL_calcUnitCuts(allVertices, unit_normal_vectors_allFaces, facesForEachVertex);

% Calculate if the vertex is convex
flag_vertexIsNonConvex_allPolys = fcn_INTERNAL_calcNonConvexVertex(facesForEachVertex_from, facesForEachVertex_to, unit_normal_vectors_allFaces);

%         % Find the unit_vertex_projection_vectors
%         if NumUniqueVerticies>2
%             pseudo_vertex_projection_vectors          = (unit_normal_vectors_thisPoly(2:end,:) + unit_normal_vectors_thisPoly(1:end-1,:))/2;
%             pseudo_vertex_projection_vectors_fullPoly = [pseudo_vertex_projection_vectors(end,:);pseudo_vertex_projection_vectors];  % vector 1 is the same as the last one
%             pseudo_vertex_projection_vector_lengths   = sum(pseudo_vertex_projection_vectors_fullPoly.^2,2).^0.5;
%             unit_vertex_projection_vectors_thisPoly   = pseudo_vertex_projection_vectors_fullPoly./pseudo_vertex_projection_vector_lengths;
%         elseif NumUniqueVerticies==2 % Line segment
%             pseudo_vertex_projection_vectors_fullPoly = unit_vectors_vertex_to_vertex_fullPoly;
%             pseudo_vertex_projection_vector_lengths   = sum(pseudo_vertex_projection_vectors_fullPoly.^2,2).^0.5;
%             unit_vertex_projection_vectors_thisPoly = pseudo_vertex_projection_vectors_fullPoly./pseudo_vertex_projection_vector_lengths;
%         elseif NumUniqueVerticies==1 % Point
%             pseudo_vertex_projection_vectors_fullPoly = zeros(size(unit_vectors_vertex_to_vertex_fullPoly));
%             unit_vertex_projection_vectors_thisPoly = pseudo_vertex_projection_vectors_fullPoly;
%         else
%             warning('on','backtrace');
%             warning('Expecting 1 or more points, but no vertices are available for vertex projection calculations! Throwing an error');
%             error('Do not know how to calculate vertex projection for one point!')
%         end
% 
%         % Find the vector direction of unit cuts
%         % See the documentation.
%         if NumUniqueVerticies>2
%             vector_sums = sum(unit_normal_vectors_thisPoly.*unit_vertex_projection_vectors_thisPoly,2);
%             d = 1./vector_sums;
%         elseif NumUniqueVerticies==2 || NumUniqueVerticies==1 % Line segment or point
%             d = 1;
%         else
%             warning('on','backtrace');
%             warning('Expecting 2 or more points, but less than 2 are available for vertex projection calculations! Throwing an error');
%             error('Do not know how to calculate vertex projection for one point!')
%         end
%         vector_direction_of_unit_cut = d.*unit_vertex_projection_vectors_thisPoly;
% 
%         % Check which ones are NOT convex. This is done by doing cross product of
%         % vectors in sequence to see if their unit normals are in same or
%         % opposite directions
%         cross_products = cross([unit_normal_vectors_thisPoly(1:end-1,:) zeros(NumUniqueVerticies,1)],[unit_normal_vectors_thisPoly(2:end,:) zeros(NumUniqueVerticies,1)],2);
%         cross_products_fullPoly = [cross_products(end,:); cross_products];
%         cross_product_results = cross_products_fullPoly(:,3);
%         flag_vertexIsNonConvex = cross_product_results<0;

% Save results
unit_normal_vectors            = unit_normal_vectors_allFaces;
vector_direction_of_unit_cut   = vector_direction_of_unit_cut_allVertices;
flag_vertexIsNonConvex         = flag_vertexIsNonConvex_allPolys;

% elseif 0==flag_useCells
%     unit_normal_vectors            = unit_normal_vectors_allFaces{1};
%     vector_direction_of_unit_cut   = vector_direction_of_unit_cut_allVertices{1};
%     flag_vertexIsNonConvex         = flag_vertexIsNonConvex_allPolys{1};
% else
%     unit_normal_vectors            = unit_normal_vectors_allFaces;
%     vector_direction_of_unit_cut   = vector_direction_of_unit_cut_allVertices;
%     flag_vertexIsNonConvex         = flag_vertexIsNonConvex_allPolys;
% 
% end

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
function [facesForEachVertex, facesForEachVertex_from, facesForEachVertex_to] = fcn_INTERNAL_findFacesForEachVertex(allVertices, allFaces)

dimension_of_points = length(allVertices(1,:));

Nvertices = length(allVertices(:,1));
LongestFace = length(allFaces(1,:));

facesForEachVertex      = nan(Nvertices, LongestFace);
facesForEachVertex_from = nan(Nvertices, 1);
facesForEachVertex_to   = nan(Nvertices, 1);

% Loop through vertices, finding which faces are involved
for ith_vertex = 1:Nvertices
    % Find all faces that use this vertex
    flag_facesWithThisVertex      = allFaces==ith_vertex;
    flag_facesWithThisVertex_from = allFaces(:,2)==ith_vertex;
    flag_facesWithThisVertex_to   = allFaces(:,1)==ith_vertex;
    
    summedFlags = sum(flag_facesWithThisVertex,2);
    facesWithThisVertex      = (find(summedFlags))';
    facesWithThisVertex_from = (find(flag_facesWithThisVertex_from))';
    facesWithThisVertex_to   = (find(flag_facesWithThisVertex_to))';
    
    if 2==dimension_of_points
        if length(facesWithThisVertex)~=2
            warning('more than 2 faces found on a vertex. This should not happen');
        end
        facesForEachVertex(ith_vertex,:)      = facesWithThisVertex;
        facesForEachVertex_from(ith_vertex,:) = facesWithThisVertex_from;
        facesForEachVertex_to(ith_vertex,:)   = facesWithThisVertex_to;
    else
        error('3D not coded yet');
    end


end % Ends for loop

end % Ends fcn_INTERNAL_findFacesForEachVertex


%% fcn_INTERNAL_calcNormalVectors
function    unit_normal_vectors_allFaces = fcn_INTERNAL_calcNormalVectors(allVertices, allFaces)


% Is this 2D or 3D?
dimension_of_points = length(allVertices(1,:));

% Calculate the unit vectors for each edge
if 2==dimension_of_points
    % Get all the indicies for all faces
    startPointIndices = allFaces(:,1);
    endPointIndicies  = allFaces(:,2);
    startVerticies = allVertices(startPointIndices,:);
    endVerticies   = allVertices(endPointIndicies,:);
   

    % find distances and unit vectors from vertex to vertex
    difference_vertex_to_vertex = endVerticies - startVerticies;
    distances_vertex_to_vertex = sum(difference_vertex_to_vertex.^2,2).^0.5;
    unit_vectors_vertex_to_vertex = difference_vertex_to_vertex./distances_vertex_to_vertex;

    % Rotate by 90 degrees to get unit normal vectors
    unit_normal_vectors_allFaces = unit_vectors_vertex_to_vertex*[0 1; -1 0];

elseif 3==dimension_of_points
    % Initialize outputs
    % Nfaces    = length(allFaces(:,1));
    % unit_normal_vectors_allFaces = nan(Nfaces, 2);


    error('Not coded yet');

    % If 3D, have to loop through faces
    % Loop through faces
    for ith_face = 1:Nfaces
        thisFace = allFaces(ith_face,:);
        thisFaceVertices = allFaces(thisFace,:);

    end
else
    warning('on','backtrace');
    warning('A vector was given that has dimension: %.0d, where 2D was expected',dimension_of_points);
    error('Function not yet coded for anything other than 2D');
end


%
%
% NumUniqueVerticies = length(thisFaceVertices(:,1))-1;
%
% % Calculate the unit vectors for each edge
% if 2==dimension_of_points
%     % find distances and unit vectors from vertex to vertex
%     difference_vertex_to_vertex = thisFaceVertices(2:end,:)-thisFaceVertices(1:end-1,:);
%     distances_vertex_to_vertex = sum(difference_vertex_to_vertex.^2,2).^0.5;
%     unit_vectors_vertex_to_vertex = difference_vertex_to_vertex./distances_vertex_to_vertex;
%
%     % Repeat last vector so it has same length as points
%     unit_vectors_vertex_to_vertex_fullPoly = [unit_vectors_vertex_to_vertex; unit_vectors_vertex_to_vertex(1,:)];
%
%     % Rotate by 90 degrees to get unit normal vectors
%     unit_normal_vectors_thisPoly = unit_vectors_vertex_to_vertex_fullPoly*[0 1; -1 0];
%
%     % Find the unit_vertex_projection_vectors
%     if NumUniqueVerticies>2
%         pseudo_vertex_projection_vectors          = (unit_normal_vectors_thisPoly(2:end,:) + unit_normal_vectors_thisPoly(1:end-1,:))/2;
%         pseudo_vertex_projection_vectors_fullPoly = [pseudo_vertex_projection_vectors(end,:);pseudo_vertex_projection_vectors];  % vector 1 is the same as the last one
%         pseudo_vertex_projection_vector_lengths   = sum(pseudo_vertex_projection_vectors_fullPoly.^2,2).^0.5;
%         unit_vertex_projection_vectors_thisPoly   = pseudo_vertex_projection_vectors_fullPoly./pseudo_vertex_projection_vector_lengths;
%     elseif NumUniqueVerticies==2 % Line segment
%         pseudo_vertex_projection_vectors_fullPoly = unit_vectors_vertex_to_vertex_fullPoly;
%         pseudo_vertex_projection_vector_lengths   = sum(pseudo_vertex_projection_vectors_fullPoly.^2,2).^0.5;
%         unit_vertex_projection_vectors_thisPoly = pseudo_vertex_projection_vectors_fullPoly./pseudo_vertex_projection_vector_lengths;
%     elseif NumUniqueVerticies==1 % Point
%         pseudo_vertex_projection_vectors_fullPoly = zeros(size(unit_vectors_vertex_to_vertex_fullPoly));
%         unit_vertex_projection_vectors_thisPoly = pseudo_vertex_projection_vectors_fullPoly;
%     else
%         warning('on','backtrace');
%         warning('Expecting 1 or more points, but no vertices are available for vertex projection calculations! Throwing an error');
%         error('Do not know how to calculate vertex projection for one point!')
%     end
%
%     % Find the vector direction of unit cuts
%     % See the documentation.
%     if NumUniqueVerticies>2
%         vector_sums = sum(unit_normal_vectors_thisPoly.*unit_vertex_projection_vectors_thisPoly,2);
%         d = 1./vector_sums;
%     elseif NumUniqueVerticies==2 || NumUniqueVerticies==1 % Line segment or point
%         d = 1;
%     else
%         warning('on','backtrace');
%         warning('Expecting 2 or more points, but less than 2 are available for vertex projection calculations! Throwing an error');
%         error('Do not know how to calculate vertex projection for one point!')
%     end
%     vector_direction_of_unit_cut = d.*unit_vertex_projection_vectors_thisPoly;
%
%     % Check which ones are NOT convex. This is done by doing cross product of
%     % vectors in sequence to see if their unit normals are in same or
%     % opposite directions
%     cross_products = cross([unit_normal_vectors_thisPoly(1:end-1,:) zeros(NumUniqueVerticies,1)],[unit_normal_vectors_thisPoly(2:end,:) zeros(NumUniqueVerticies,1)],2);
%     cross_products_fullPoly = [cross_products(end,:); cross_products];
%     cross_product_results = cross_products_fullPoly(:,3);
%     flag_vertexIsNonConvex = cross_product_results<0;
% else
%     warning('on','backtrace');
%     warning('A vector was given that has dimension: %.0d, where 2D was expected',dimension_of_points);
%     error('Function not yet coded for anything other than 2D');
% end

end % Ends fcn_INTERNAL_calcNormalVectors

%% fcn_INTERNAL_calcUnitCuts
function vector_direction_of_unit_cut_allVertices = fcn_INTERNAL_calcUnitCuts(allVertices, unit_normal_vectors_allFaces, facesForEachVertex)

% Is this 2D or 3D?
dimension_of_points = length(allVertices(1,:));

% Initialize outputs
vector_direction_of_unit_cut_allVertices = nan(length(allVertices(:,1)),dimension_of_points);

% Loop through vertices, finding which faces are involved
for ith_vertex = 1:length(allVertices)
    
    facesWithThisVertex = (facesForEachVertex(ith_vertex,:))';
    
    if 2==dimension_of_points

        % Set up A matrix for linear equation solution
        Amatrix = unit_normal_vectors_allFaces(facesWithThisVertex,:);

        if any(isnan(Amatrix),'all') || any(isinf(Amatrix),'all') || rank(Amatrix)<dimension_of_points
            % error('Degenerate A matrix found');
            % Do nothing - returns NaN by default 
        else
            % Solve linear equation
            vector_direction_of_unit_cut_allVertices(ith_vertex,:) = (Amatrix\ones(dimension_of_points,1))';
        end
    else
        error('3D not coded yet');
    end


end % Ends for loop

end % Ends fcn_INTERNAL_calcUnitCuts

%% fcn_INTERNAL_calcNonConvexVertex
function flag_vertexIsNonConvex = fcn_INTERNAL_calcNonConvexVertex(facesForEachVertex_from, facesForEachVertex_to, unit_normal_vectors_allFaces)

% Is this 2D or 3D?
dimension_of_points = length(unit_normal_vectors_allFaces(1,:));
Nvertices = length(facesForEachVertex_from(:,1));

if 2==dimension_of_points
    % Find all faces for this vertex
    startFaces = facesForEachVertex_from(:,1);
    endFaces   = facesForEachVertex_to(:,1);

    % Find their cross product
    cross_products = cross([unit_normal_vectors_allFaces(startFaces,:) zeros(Nvertices,1)],[unit_normal_vectors_allFaces(endFaces,:) zeros(Nvertices,1)],2);
    cross_product_results = cross_products(:,3);

    % Vertex is non convex if cross product less than zero
    flag_vertexIsNonConvex = cross_product_results<0;

else
    error('3D not coded yet');
end



end % Ends fcn_INTERNAL_calcNonConvexVertex