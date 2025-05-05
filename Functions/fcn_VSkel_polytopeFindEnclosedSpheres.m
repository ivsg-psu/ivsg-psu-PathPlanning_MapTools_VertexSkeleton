function [sphereEdgeRadii, definingEdges] = ...
    fcn_VSkel_polytopeFindEnclosedSpheres(vertices, unit_normal_vectors, unit_vertex_projection_vectors, vector_direction_of_unit_cut, flag_vertexIsNonConvex, max_edge_cuts, varargin)

%% fcn_VSkel_polytopeFindEnclosedSpheres
% finds the enclosed sphere radii for each vertex point, to each other edge
%
% FORMAT:
%
% [sphereRadii, definingEdges] = ...
%     fcn_VSkel_polytopeFindEnclosedSpheres(vertices, unit_normal_vectors, unit_vertex_projection_vectors, flag_vertexIsNonConvex, (fig_num))
%
% INPUTS:
%
%     vertices: a (M+1)-by-2 matrix of xy points with row1 = rowm+1, where
%         M is the number of the individual polytope vertices
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
%     sphereEdgeRadii: a cell array of dimension N containing, in each
%     cell, an array of radii that cause that vertex projection to contact
%     an edge. In each cell array, the are Mx1 radii vectors, where M is =
%     N-2 and N is the number of vertices. The radii are ordered so that
%     the first radii cell array corresponds to the first vertex, etc.
%
%     definingEdges: a cell array of dimension N containing, in
%     each cell, an array of which edges define each vertex. In each cell
%     array, there are Mx1 defining edges, where M is = N-2 and N is the
%     number of vertices. The defining edges match the radii ordering, e.g.
%     vertex 2's 3rd sphereRadii edge interaction ID will be in cell array
%     2, in the 3rd row.
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

% Is this 2D or 3D?
dimension_of_points = length(vertices(1,:));
NumUniqueVerticies = length(vertices(:,1))-1;

% Calculate the unit vectors for each edge
if 2==dimension_of_points
    % Initialize variables
    sphereEdgeRadii      = cell(NumUniqueVerticies+1,1);
    sphereEdgeCenters    = cell(NumUniqueVerticies+1,1);
    definingEdges        = cell(NumUniqueVerticies+1,1);

    % For each vertex, solve for the radii to all the non-participating
    % edges. Non-participating edges are those that are not to either side
    % of this vertex.
    for ith_vertex = 1:NumUniqueVerticies

        % Find which edges are intersecting with the current vertex        
        [radiiFromVertexToEdge, sphereEdgeCenterArray, edgesConstrainingRadii] = ...
            fcn_INTERNAL_checkVertexProjectionsOntoEdges(vertices, ith_vertex, unit_normal_vectors, unit_vertex_projection_vectors, vector_direction_of_unit_cut, max_edge_cuts);

        % % % Find which verticies are intersecting with the current vertex        
        %  [radiiFromVertexToVertex, sphereVertexCenterArray, verticesConstrainingRadii] = ...
        %      fcn_INTERNAL_checkVertexProjectionsOntoVertices(vertices, ith_vertex, unit_normal_vectors, vertexProjection, flag_vertexIsNonConvex);


        % Save matricies to cell arrays for this vertex
        sphereEdgeRadii{ith_vertex}     = radiiFromVertexToEdge;
        sphereEdgeCenters{ith_vertex}   = sphereEdgeCenterArray;
        definingEdges{ith_vertex}       = edgesConstrainingRadii;
        % sphereVertexRadii{ith_vertex}   = radiiFromVertexToVertex;
        % sphereVertexCenters{ith_vertex} = sphereVertexCenterArray;
        % definingVerticies{ith_vertex}   = verticesConstrainingRadii;
    end

    % Repeat last value as first, since this is the repeated 1st vertex
    sphereEdgeRadii{end}     = sphereEdgeRadii{1};
    sphereEdgeCenters{end}   = sphereEdgeCenters{1};
    definingEdges{end}       = definingEdges{1};

else
    warning('on','backtrace');
    warning('A vector was given that has dimension: %.0d, where 2D was expected',dimension_of_points);
    error('Function not yet coded for anything other than 2D');
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
    
    tiledlayout('flow');

    % Find size of vertex domain
    max_XY = max(vertices);
    min_XY = min(vertices);
    sizePlot = max(max_XY) - min(min_XY);
    nudge = sizePlot*0.006;


    for this_vertex = 1:NumUniqueVerticies

        radiiFromVertexToEdge  = sphereEdgeRadii{this_vertex};
        sphereEdgeCenterArray  = sphereEdgeCenters{this_vertex};
        edgesConstrainingRadii = definingEdges{this_vertex};
     
        nexttile;

        fcn_VSkel_plotPolytopeDetails(...
            vertices,...
            (unit_normal_vectors), ...  % unit_normal_vectors
            ([]), ...  % unit_vertex_projection_vectors
            ([]), ... % vector_direction_of_unit_cut
            (flag_vertexIsNonConvex),...  % flag_vertexIsNonConvex
            (1),...  % flag_plotEdgeGhostlines
            (1),...  % flag_plotVertexProjectionGhostlines
            (fig_num));  % fig_num

        % Keep the axis from this standard plot, so that all subplots are
        % same
        goodAxis = axis;

        % Plot the current vertex with a red circle
        plot(vertices(this_vertex,1), vertices(this_vertex,2),'ro','MarkerSize',3);

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
                sprintf('%.0dto %.0d',this_vertex,circleEdge), 'Color',colorUsed);
        end

        % Make axis back to before
        axis(goodAxis);
        title(sprintf('Vertex: %.0d',this_vertex));
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
methodToCheck = 1;
tolerance = 1E-5;

edgeLength     = real(sum((edgeEnd - edgeStart).^2,2).^0.5);

% For debugging

if 1==flag_do_debug
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
    if dotProduct<(0-tolerance) || dotProduct>(edgeLength+tolerance)
        isOnEdge = 0;
    else
        % Do dot product to see if not orthogonal
        ortho_edgeVector = unit_edgeVector*[0 1; -1 0];
        dotProduct = sum(ortho_edgeVector.*testVector,2);
        if dotProduct<(0-tolerance) || dotProduct>(0+tolerance)
            isOnEdge = 0;
        else
            isOnEdge = 1;
        end
    end

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

if 1==flag_do_debug
    figure(debug_fig_num);
    if 1==isOnEdge
        plot(testPoint(:,1),testPoint(:,2),'o','MarkerSize',20,'Color',colorUsed);
    else
        plot(testPoint(:,1),testPoint(:,2),'x','MarkerSize',20,'Color',[1 0 0]);
    end
end

end % Ends fcn_INTERNAL_isPointOnEdge

