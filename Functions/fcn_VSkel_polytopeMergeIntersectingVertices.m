function vertexSkeletonStartingPolytopes = ...
    fcn_VSkel_polytopeMergeIntersectingVertices( ...
    vertices, ...
    flag_vertexIsNonConvex,...
    bounaryEngagedAtMinCut, ...
    indices_repeated, ...
    intersection_points, ...
    flag_vertexIsNonConvex, ...
    varargin) % definingBoundaries, unit_normal_vectors, unit_vertex_projection_vectors, vector_direction_of_unit_cut, 

%% fcn_VSkel_polytopeMergeIntersectingVertices
% finds the minimum cut possible before a vertex is eliminated
%
% FORMAT:
%
% [min_cut, indices_repeated] = ...
%     fcn_VSkel_polytopeMergeIntersectingVertices(vertices, ...
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
%     vertexSkeletonStartingPolytopes: a structure defining the vertex
%     starting points for the next vertex level, e.g. the following fields
%     are filled in for each polytope:
%
%     vertexSkeletonStartingPolytopes.polytope(polytope_index).vertices
%     
%     A vertex level has the following structure:
%
%        vertexSkeleton(depth_index).FIELDS
%
%     where depth_index refers to the index of the depth of cut. The
%     subfields are indexed by the polytope_index, which counts the number
%     of active polytopes at this given cut depth:
%
%        vertexSkeleton(depth_index).polytope(polytope_index).SUBFIELDS
% 
%     where subfields are as follows:
%
%        vertexSkeleton(depth_index).polytope(polytope_index).vertices                       = [1 2; 3 4; 5 6; 1 2]; % for 2D, the 1st and last points are same, so this is of size N+1 if there are N unique vertices 
%        vertexSkeleton(depth_index).polytope(polytope_index).unit_vertex_projection_vectors = [1 2; 3 4; 5 6];      % same number of rows as vertices; 
%        vertexSkeleton(depth_index).polytope(polytope_index).vector_direction_of_unit_cut = [1 2; 3 4; 5 6]         % same number of rows as vertices;
%        vertexSkeleton(depth_index).polytope(polytope_index).flag_vertexIsNonConvex = [0; 0; 0; 0];                 % same number of rows as vertices;
%        vertexSkeleton(depth_index).polytope(polytope_index).boundaries{boundary_index}   = [1 2; 3 4];             % number of bounding edges. In 3D, the number of rows can be large but must start/end on same last point 
%        vertexSkeleton(depth_index).polytope(polytope_index).unit_normal_vector{edge_index} = [1 2];                % normal vector to the boundary, so will be a [1xD] where D is the dimension;
%
% DEPENDENCIES:
%
%     fcn_DebugTools_checkInputsToFunctions
%
% % EXAMPLES:
% For additional examples, see: script_test_fcn_VSkel_polytopeMergeIntersectingVertices
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

        URHERE - fix argument list, docstrings, and checks

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

% Initialize which vertices have been eliminated
indices_repeated_alreadyAnalyzed = indices_repeated*0;
% Calculate the unit vectors for each edge
if 2==dimension_of_points

    %%%%%%
    % Steps:
    % 1) Loop through all indices_repeated, making sure not to analyze
    % indicies previously treated, check if intersection points are all
    % convex. If yes, then flag that vertices are all merged into one. Save
    % result as 
    % 
    %    flag_verticesToMerge
    % 
    % with 0 meaning no mergers, 1 meaning the first merger, 2 meaning the
    % second merger, 
    % etc. For example, flag_verticesToMerge would be, if vertex 4 is the
    % first merger and is merging with 5 and 6 but not 7. [0 0 0 1 1 1 0 0 ...]
    %
    % 2) Making sure not to analyze indicies previously treated, remaining
    % indicies must not all be convex. This indicates that verticies are
    % separating, and new polytopes must be made at each of these vertices.
    % This means that verticies have to be grouped into sub-polytopes. To
    % do this, flag all "break" points using
    % 
    %    flag_verticiesWherePolytopesSeparate
    %
    % with 0 meaning no break points, 1
    % meaning break point 1, 2 for break point 2, etc. The number of new
    % polytopes to create will be the max. As well, update the list of
    % point indicies and insert a -N into the list at the location of break
    % point, where N indicates the location where vertex N intrudes. So the
    % updatedVertexSequence will have values such as [1; 2; 3; -7; 4; 5;
    % etc.]. This would indicate that vertex 7 intrudes between vertices 3
    % and 4 in the resulting polygon.
    %
    % 3) Create a cell array to store which vertices will go into which
    % polytopes. Starting from the first vertex, proceed along and check
    % for mergers. If there is a merger at this vertex and no separations,
    % drop the points other than the lowest index. If there is a separation
    % at this index, use the separation point to close off the previous
    % polygon and start a new one.

    % Initialize output
    vertexSkeletonStartingPolytopes = struct;

    %%%%%%
    % STEP 1) Loop through all indices_repeated, making sure not to analyze
    % indicies previously treated, check if intersection points are all
    % convex. If yes, then flag that vertices are all merged into one. Save
    % result as 
    % 
    %    flag_verticesToMerge
    % 
    % with 0 meaning no mergers, 1 meaning the first merger, 2 meaning the
    % second merger, 
    % etc. For example, flag_verticesToMerge would be, if vertex 4 is the
    % first merger and is merging with 5 and 6 but not 7. [0 0 0 1 1 1 0 0 ...]

    flag_verticesToMerge = zeros(NumUniqueVerticies,1);
    Num_mergers = 0;
    mergedVertexIDList = (1:NumUniqueVerticies)';
    for ith_repeatTest = 1:length(indices_repeated)
        % Make sure not to analyze indices that were previously analyzed
        if indices_repeated_alreadyAnalyzed(ith_repeatTest)==0
            this_index = indices_repeated(ith_repeatTest);

            indices_in_this_merge = fcn_INTERNAL_findIndicesInThisIntersection(vertices, indices_repeated, this_index);
           
            % Set indicies in for loop ahead of this point to not repeat
            % the search
            indices_repeated_alreadyAnalyzed(indices_in_this_merge) = 1;
            
            % Check if these vertices are all convex
            flag_thisMergeIsAllConvex = all(flag_vertexIsNonConvex(indices_in_this_merge)==0);

            if flag_thisMergeIsAllConvex
                Num_mergers = Num_mergers+1;
                flag_verticesToMerge(indices_in_this_merge,:) = Num_mergers;
                mergedVertexIDList(indices_in_this_merge) = min(indices_in_this_merge);
            end
        end % Ends test to see if this is a repeated analysis
    end % Ends for loop through repeated indices

    %%%%%%%%%%%%%%%%%%%%
    % STEP 2) Making sure not to analyze indicies previously treated, remaining
    % indicies must not all be convex. This indicates that verticies are
    % separating, and new polytopes must be made at each of these vertices.
    % This means that verticies have to be grouped into sub-polytopes. We
    % first increment the number of separations. As well, update the list
    % of point indicies and insert a -N into the list at the location of
    % break point, where N indicates the location where vertex N intrudes.
    % So the updatedVertexSequence will have values such as [1; 2; 3; -7;
    % 4; 5; etc.]. This would indicate that vertex 7 intrudes between
    % vertices 3 and 4 in the resulting polygon.

    flag_verticiesWherePolytopesSeparate = [];
    Num_separations = 0;
    previousIndex = 1;
    for ith_repeatTest = 1:length(indices_repeated)
        % Make sure not to analyze indices that were previously analyzed
        if indices_repeated_alreadyAnalyzed(ith_repeatTest)==0
            this_index = indices_repeated(ith_repeatTest);
            % Update the sequence of points
            flag_verticiesWherePolytopesSeparate = [flag_verticiesWherePolytopesSeparate; mergedVertexIDList(previousIndex:this_index,1)]; %#ok<AGROW>
            previousIndex = this_index;

            indices_in_this_merge = fcn_INTERNAL_findIndicesInThisIntersection(vertices, indices_repeated, this_index);

            % Set indicies in for loop ahead of this point to not repeat
            % the search
            indices_repeated_alreadyAnalyzed(indices_in_this_merge) = 1;
            
            % Confirm at least one of thse is NOT convex
            flag_thisMergeIsNOTAllConvex = any(flag_vertexIsNonConvex(indices_in_this_merge)==1);

            if flag_thisMergeIsNOTAllConvex
                Num_separations = Num_separations+1;
                flag_verticesToMerge(indices_in_this_merge,:) = Num_mergers;
            else
                warning('on','backtrace');
                warning('Indicies were tagged as non-convex but were convex. Not sure how to proceed?');
                error('Convex indices found when non-convex indices expected? Exiting.');
            end
        end % Ends test to see if this is a repeated analysis
    end % Ends for loop through repeated indices

    % Update the sequence of points
    if this_index~=NumUniqueVerticies
        flag_verticiesWherePolytopesSeparate = [flag_verticiesWherePolytopesSeparate; mergedVertexIDList(prevous_index:NumUniqueVerticies,1)];
    end

    disp(flag_verticiesWherePolytopesSeparate);


    %%%%%
    % Step 3: Create a cell array to store which vertices will go into which
    % polytopes. Starting from the first vertex, proceed along and check
    % for mergers. If there is a merger at this vertex and no separations,
    % drop the points other than the lowest index. If there is a separation
    % at this index, use the separation point to close off the previous
    % polygon and start a new one.
    for ith_polytope = 1:(Num_separations+1)

        vertexSkeletonStartingPolytopes.polytope(polytope_index).vertices
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
        (fig_num));  % fig_num

    % Plot the intersection_points
    plot(intersection_points(:,1), intersection_points(:,2),'.','MarkerSize',30);


    % Plot the spheres, and label their centers
    for this_repeat = 1:length(indices_repeated)
        ith_vertex = indices_repeated(this_repeat);
        circleCenter = vertices(ith_vertex,:) + vector_direction_of_unit_cut(ith_vertex,:)*min_cut;
        circleRadius = min_cut;
        circleEdge = bounaryEngagedAtMinCut(ith_vertex,1);

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

%% fcn_INTERNAL_findIndicesInThisIntersection
function indices_in_this_merge = fcn_INTERNAL_findIndicesInThisIntersection(vertices, indices_repeated, this_index)
% Given a list of vertices, a list of possible repeats, and this index,
% checks to see if the repeated points are the SAME point, within a
% tolerance. For those that are the same, returns the indices that are all
% the same.

% Note the possible repeats
possibleRepeatVertices = vertices(indices_repeated,:);

% Check if any other intersection points are the same
% NOTE: may need to add a tolerance here
current_point = vertices(this_index,:);
differences_with_current = sum((possibleRepeatVertices - current_point).^2,2);
tolerance = 1E-6;
flag_vertices_same = (differences_with_current<tolerance);

% Find which vertices are intersecting here
indices_in_this_merge = indices_repeated(flag_vertices_same);

end % Ends fcn_INTERNAL_findIndicesInThisIntersection