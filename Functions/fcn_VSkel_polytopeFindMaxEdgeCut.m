function max_edge_cuts = fcn_VSkel_polytopeFindMaxEdgeCut(vertices, unit_normal_vectors, unit_vertex_projection_vectors, varargin)

%% fcn_VSkel_polytopeFindMaxEdgeCut
% for each edge in a 2D polytope, finds the maximum cut that is possible
% before the edge disappears. This is done by examining when the vertex
% projection vectors on each side of an edge converge to each other.
% Returns Inf values if they do not intersect.
%
% FORMAT:
%
% max_edge_cuts = fcn_VSkel_polytopeFindMaxEdgeCut(vertices, unit_vertex_projection_vectors, (fig_num))
%
% INPUTS:
%
%     vertices: a (M+1)-by-2 matrix of xy points with row1 = rowm+1, where
%         M is the number of the individual polytope vertices
%
%     unit_normal_vectors: an (M+1)-by-2 matrix of the unit vectors that
%     point inward as measured from one vertex to the next. The vector is
%     assumed to be attached to the start of the edge given by vertex M.
%
%     unit_vertex_projection_vectors: an (M+1)-by-2 matrix of the unit
%     vectors that project in the direction that each of the M verticies
%     will move. The vector is assumed to be attached to the start of the
%     edge given by vertex M.
%
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
%     max_edge_cuts: a (M+1)-by-1 matrix of distances, for each vertex,
%     that can be cut from the edges before the edge "self-intersects",
%     e.g. the projection from each side of the edges meet. If there is not
%     any intersection, the distance is infinite.
%
%
% DEPENDENCIES:
%
%     fcn_DebugTools_checkInputsToFunctions
%
% % EXAMPLES:
% For additional examples, see: script_test_fcn_VSkel_polytopeFindMaxEdgeCut
%
% This function was written on 2025_05_04 by S. Brennan
% Questions or comments? sbrennan@psu.edu
%

% Revision History:
% 2025_05_04 - S. Brennan
% -- first write of code

% TO DO
% -- none

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==4 && isequal(varargin{end},-1))
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
        narginchk(3,4);

        % Check the vertices input
        fcn_DebugTools_checkInputsToFunctions(...
            vertices, '2or3column_of_numbers');

        NumUniqueVerticies = length(vertices(:,1));
        
        % Check the unit_normal_vectors input
        fcn_DebugTools_checkInputsToFunctions(...
            unit_normal_vectors, '2or3column_of_numbers',NumUniqueVerticies);
        
        % Check the vertex_projection_vectors input
        fcn_DebugTools_checkInputsToFunctions(...
            unit_vertex_projection_vectors, '2or3column_of_numbers',NumUniqueVerticies);
    end
end

% Does user want to show the plots?
flag_do_plot = 0; % Default is no plotting
if  4 == nargin && (0==flag_max_speed) % Only create a figure if NOT maximizing speed
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

% Is this 2D or 3D
dimension_of_points = length(vertices(1,:));
NumUniqueVerticies = length(vertices(:,1))-1;



% Find start/end vectors for each edge
edgeStartVectors = unit_vertex_projection_vectors(1:NumUniqueVerticies,:);
edgeEndVectors   = unit_vertex_projection_vectors(2:NumUniqueVerticies+1,:);

% Initialize outputs
max_edge_cuts = inf(NumUniqueVerticies+1,1);

% Calculate the unit vectors for each edge
if 2==dimension_of_points

    % Find size of vertex domain
    max_XY = max(vertices);
    min_XY = min(vertices);
    sizePlot = 2*(max(max_XY) - min(min_XY));
    
    % Use cross product to find all edges that have finite intersection
    % points
    cross_result = cross([edgeStartVectors zeros(NumUniqueVerticies,1)],[edgeEndVectors zeros(NumUniqueVerticies,1)],2);
    flag_edgesHaveIntersections = cross_result(:,3)>0;
    edgesToCheck = find(flag_edgesHaveIntersections);

    for ith_search = 1:length(edgesToCheck)
        thisEdgeID = edgesToCheck(ith_search);
        vertexStart = vertices(thisEdgeID,:);
        vertexEnd   = vertices(thisEdgeID+1,:);
        projectionStart = sizePlot*edgeStartVectors(thisEdgeID,:);
        projectionEnd   = sizePlot*edgeEndVectors(thisEdgeID,:);
        edgeUnitNormal = unit_normal_vectors(thisEdgeID,:);
 

        % For debugging
        if 1==0
            fcn_VSkel_plotPolytopeDetails(...
                [vertexStart; vertexEnd; vertexStart],...
                ([]), ...  % unit_normal_vectors
                ([projectionStart; projectionEnd; projectionStart]), ...
                ([]), ... % vector_direction_of_unit_cut
                ([]),...  % flag_vertexIsNonConvex
                (8888));
        end

        % Call the path library to find intersections
        % FORMAT:
        % [distance, location, path_segment, t, u] = ...
        %     fcn_Path_findProjectionHitOntoPath(path,...
        %     sensor_vector_start,sensor_vector_end,...
        %     (flag_search_type),(fig_num))
        path = [vertexStart; vertexStart+projectionStart];
        [distance, location] = ...
            fcn_Path_findProjectionHitOntoPath(path,...
            vertexEnd,vertexEnd+projectionEnd,...
            (0),([]));
        if ~isnan(distance)
            projectionVector = location - vertexStart;
            % Do the dot product to find the cutDistance
            cutDistance = sum((edgeUnitNormal.*projectionVector),2);
            max_edge_cuts(thisEdgeID) = cutDistance;
        end          
    end
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
    tempPlotLength = max_edge_cuts;
    tempPlotLength(isinf(max_edge_cuts)) = 0;
    valuesToPlot = unit_normal_vectors.*tempPlotLength;

    fcn_VSkel_plotPolytopeDetails(...
       vertices,...
       (valuesToPlot), ...  % unit_normal_vectors
       (unit_vertex_projection_vectors), ...  % unit_vertex_projection_vectors
       (unit_vertex_projection_vectors), ... % vector_direction_of_unit_cut
       ([]),...  % flag_vertexIsNonConvex
       (1),...  % flag_plotEdgeGhostlines
       (1),...  % flag_plotVertexProjectionGhostlines
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