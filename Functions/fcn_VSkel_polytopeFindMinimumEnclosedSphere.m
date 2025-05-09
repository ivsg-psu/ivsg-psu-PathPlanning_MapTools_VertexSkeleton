function [min_cut, boundaryEngagedAtMinCut, ...
    indices_repeated, intersection_points] = ...
    fcn_VSkel_polytopeFindMinimumEnclosedSphere(vertices, ...
    sphereRadii, definingBoundaries, unit_normal_vectors, unit_vertex_projection_vectors, ...
    vector_direction_of_unit_cut, varargin)

%% fcn_VSkel_polytopeFindMinimumEnclosedSphere
% finds the minimum cut possible before a vertex is eliminated
%
% FORMAT:
%
% [min_cut, indices_repeated] = ...
%     fcn_VSkel_polytopeFindMinimumEnclosedSphere(vertices, ...
%     sphereRadii, definingBoundaries, unit_normal_vectors, unit_vertex_projection_vectors, ...
%     vector_direction_of_unit_cut, (fig_num)
%
% INPUTS:
%
%     vertices: a (M+1)-by-2 matrix of xy points with row1 = rowm+1, where
%         M is the number of the individual polytope vertices
%
%     sphereRadii: a cell array of dimension N containing, in each cell, an
%     array of radii. In each cell array, the are Mx1 radii vectors, where
%     M is = N-2 and N is the number of vertices. The radii are ordered so
%     that the first radii cell array corresponds to the first vertex, etc.
%
%     definingBoundaries: a cell array of dimension N containing, in each cell,
%     an array of which edges constrain each radius of each vertex. In each
%     cell array, there are Mx1 defining edges, where M is = N-2 and N is
%     the number of vertices. The defining edges match the radii ordering,
%     e.g. vertex 2's 3rd sphereRadii edge interaction ID will be in cell
%     array 2, in the 3rd row.
%
%     unit_normal_vectors: a cell array of dimension M, where
%     each index 1:M stores a N x 2 array of the unit vectors that point
%     inward as measured from one vertex to the next.
%
%     unit_vertex_projection_vectors: a cell array of M, where each index
%     1:M stores a N x 2 array of the unit vectors that point away from the
%     vertices into the nested shape inside, with M = 1 being the starting
%     unit vectors and N being smaller and smaller for each M value.
%
%     vector_direction_of_unit_cut: a cell array of dimension M, where each
%     index 1:M stores a N x 2 array of the vectors that define the
%     magnitude and diretion of the vertices movement into the nested shape
%     inside, assuming a unit magnitude cut. 
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
%     min_cut: the smallest cut possible among all the edges
%
%     boundaryEngagedAtMinCut: a (M+1)-by-2 matrix of the ID of the boundary
%     constraining the radius of each vertex.
%
%     indices_repeated: which vertex indicies have repeated cut depths
%     that match the minimum
%
%     intersection_points: locations of all verticies after the min_cut is
%     applied
%
% DEPENDENCIES:
%
%     fcn_DebugTools_checkInputsToFunctions
%
% % EXAMPLES:
% For additional examples, see: script_test_fcn_VSkel_polytopeFindMinimumEnclosedSphere
%
% This function was written on 2025_05_03 by S. Brennan
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

% To-DO
% 2025-05-05 - need to remove unused arguments

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

        % Check the vertices input
        fcn_DebugTools_checkInputsToFunctions(...
            vertices, '2or3column_of_numbers');

        % Check the unit_normal_vectors input
        fcn_DebugTools_checkInputsToFunctions(...
            unit_normal_vectors, '2or3column_of_numbers',NumUniqueVerticies);

        % Check the unit_vertex_projection_vectors input
        fcn_DebugTools_checkInputsToFunctions(...
            unit_vertex_projection_vectors, '2or3column_of_numbers',NumUniqueVerticies);

        % Check the vector_direction_of_unit_cut input
        fcn_DebugTools_checkInputsToFunctions(...
            vector_direction_of_unit_cut, '2or3column_of_numbers',NumUniqueVerticies);

        assert(length(sphereRadii)==NumUniqueVerticies);
        assert(length(definingBoundaries)==NumUniqueVerticies);

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
    minCutsEachVertex   = zeros(NumUniqueVerticies,1);
    boundariesEngagedAtMinCut = zeros(NumUniqueVerticies+1,1);
    % For each vertex, solve for the radii to all the non-participating
    % edges. Non-participating edges are those that are not to either side
    % of this vertex.
    for ith_vertex = 1:NumUniqueVerticies
        this_vertex_radii = sphereRadii{ith_vertex};
        this_vertex_edges = definingBoundaries{ith_vertex};
        [this_vertex_min_cut, this_vertex_min_cut_index]   = min(this_vertex_radii);
        minCutsEachVertex(ith_vertex,1) = this_vertex_min_cut;
        boundariesEngagedAtMinCut(ith_vertex,1) = this_vertex_edges(this_vertex_min_cut_index,1);
    end
    boundariesEngagedAtMinCut(NumUniqueVerticies+1,1) = boundariesEngagedAtMinCut(1,1);

    % Find the minimum cut
    [min_cut, ~] = min(minCutsEachVertex);

    % Find repeats that have the same cut distance
    indices_repeated = find(minCutsEachVertex<(min_cut+1E5*eps));

    % Find intersection points
    intersection_points = vertices + vector_direction_of_unit_cut*min_cut;

    % Set boundaryEngagedAtMinCut
    boundaryEngagedAtMinCut = boundariesEngagedAtMinCut(indices_repeated);

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

    % check whether the figure already has data
    temp_h = figure(fig_num);
    flag_rescale_axis = 0;
    if isempty(get(temp_h,'Children'))
        flag_rescale_axis = 1;
    end

    tiledlayout('flow');

    % Find size of vertex domain
    max_XY = max(vertices);
    min_XY = min(vertices);
    sizePlot = max(max_XY) - min(min_XY);
    nudge = sizePlot*0.006;

    if flag_rescale_axis
        axis_range_x = max_XY(1,1)-min_XY(1,1);
        axis_range_y = max_XY(1,2)-min_XY(1,2);
        percent_larger = 0.3;
        axis([min_XY(1,1)-percent_larger*axis_range_x, max_XY(1,1)+percent_larger*axis_range_x,  min_XY(1,2)-percent_larger*axis_range_y, max_XY(1,2)+percent_larger*axis_range_y]);
        goodAxis = axis;
    end


    fcn_VSkel_plotPolytopeDetails(...
        vertices,...
        (unit_normal_vectors), ...  % unit_normal_vectors
        (unit_vertex_projection_vectors), ...  % unit_vertex_projection_vectors
        ([]), ... % vector_direction_of_unit_cut
        ([]),...  % flag_vertexIsNonConvex
        (1),...  % flag_plotEdgeGhostlines
        (1),...  % flag_plotVertexProjectionGhostlines
        ([]),...  % plot_formatting
        (fig_num));  % fig_num

    % Plot the intersection_points
    plot(intersection_points(:,1), intersection_points(:,2),'.','MarkerSize',30);


    % Plot the spheres, and label their centers
    for this_repeat = 1:length(indices_repeated)
        ith_vertex = indices_repeated(this_repeat);
        circleCenter = vertices(ith_vertex,:) + vector_direction_of_unit_cut(ith_vertex,:)*min_cut;
        circleRadius = min_cut;
        circleEdge = boundariesEngagedAtMinCut(ith_vertex,1);

        % Plot the circle center as a large dot, and store the color
        h_fig = plot(circleCenter(1,1),circleCenter(1,2),'.','MarkerSize',20);
        colorUsed = get(h_fig,'Color');

        % Plot the circle boundary in same color
        fcn_geometry_plotCircle(circleCenter,circleRadius,colorUsed,fig_num);

        text(circleCenter(1,1)+nudge, circleCenter(1,2)+4*nudge*(this_repeat-1),...
            sprintf('%.0dto %.0d',ith_vertex,circleEdge), 'Color',colorUsed);
    end



    % Make axis slightly larger?
    if flag_rescale_axis
        axis(goodAxis);
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
