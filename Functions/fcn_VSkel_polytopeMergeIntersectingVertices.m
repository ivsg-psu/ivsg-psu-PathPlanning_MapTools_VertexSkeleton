function vertexSkeletonStartingPolytopes = ...
    fcn_VSkel_polytopeMergeIntersectingVertices( ...
    vertices, ...
    flag_vertexIsNonConvex,...
    boundaryEngagedAtMinCut, ...
    indices_repeated, ...
    intersection_points, ...
    varargin)

%% fcn_VSkel_polytopeMergeIntersectingVertices
% separates vertices into different polygons based on intersections of
% points during cutting.
%
% FORMAT:
%
% vertexSkeletonStartingPolytopes = ...
%     fcn_VSkel_polytopeMergeIntersectingVertices( ...
%     vertices, ...
%     flag_vertexIsNonConvex,...
%     boundaryEngagedAtMinCut, ...
%     indices_repeated, ...
%     intersection_points, ...
%     (fig_num));
%
% INPUTS:
%
%     vertices: a (M+1)-by-2 matrix of xy points with row1 = rowm+1, where
%         M is the number of the individual polytope vertices
%
%     flag_vertexIsNonConvex: an N x 1 array of flags (true or false) that
%     indicate whether the vertex is not convex (1 = NOT convex)
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
% 2025_05_03 - S. Brennan
% -- first write of code

% To-DO
% 2025-05-05 
% - need to remove unused arguments (vertices)
% - merge the steps 1 and 2 to use same for loop

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==6 && isequal(varargin{end},-1))
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
        narginchk(5,6);

        % Check the vertices input
        fcn_DebugTools_checkInputsToFunctions(...
            vertices, '2or3column_of_numbers');

        NumVertices = length(vertices(:,1));

        % Check the flag_vertexIsNonConvex input
        fcn_DebugTools_checkInputsToFunctions(...
            flag_vertexIsNonConvex*1.00, '1column_of_numbers', NumVertices);

        % Check the boundaryEngagedAtMinCut input
        fcn_DebugTools_checkInputsToFunctions(...
            boundaryEngagedAtMinCut, '1column_of_numbers');

        NumIndicesRepeated = length(boundaryEngagedAtMinCut(:,1));

        % Check the indices_repeated input
        fcn_DebugTools_checkInputsToFunctions(...
            indices_repeated, '1column_of_numbers', NumIndicesRepeated);

        % Check the intersection_points input
        fcn_DebugTools_checkInputsToFunctions(...
            intersection_points, '2or3column_of_numbers',NumVertices);

    end
end

% Does user want to show the plots?
flag_do_plot = 0; % Default is no plotting
if  6 == nargin && (0==flag_max_speed) % Only create a figure if NOT maximizing speed
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
dimension_of_points = length(intersection_points(1,:));

% Calculate the unit vectors for each edge
if 2==dimension_of_points

    %%%%%%
    % Steps:
    % STEP 1) Loop through all indices_repeated, making sure not to analyze
    % indicies previously treated, find all points in this intersection.
    % Check if all these intersection points are all same. If so, flag
    % that vertices are all merged into one. Save result as
    % 
    %    flag_verticesToMerge
    % 
    % with 0 meaning no mergers, 1 meaning the first merger, 2 meaning the
    % second merger, 
    % etc. For example, flag_verticesToMerge would be, if vertex 4 is the
    % first merger and is merging with 5 and 6 but not 7, 8, 9, etc. then
    % flag_verticesToMerge would be [0 0 0 1 1 1 0 0 ...]. Finish by
    % flagging the merged points as treated.
    
    % STEP 2) Of indicies_repeated, making sure not to analyze non-convex
    % indicies previously treated, loop through each one, finding insertion
    % points for each indicating where polytopes are separating. New polytopes
    % must be made at each of these vertices. This means that verticies have to
    % be grouped into sub-polytopes in a later step. This is difficult, because
    % the grouping of points into polytopes means all indicies - even ones not
    % yet processed - need to be accounted for, otherwise polytopes can be
    % created that would then have to themselves be broken again in later
    % steps. Thus, the point of this step is to process the insertion points of
    % all vertices, thus defining where polytope "breaks" occur from one
    % polytope to another.  This is done by creating and updating a list of
    % point indicies and insert a -N into the list at the location of break
    % point, where N indicates the location where nonConvex vertex N intrudes.
    % So the updatedVertexSequence will have values such as [1; 2; 3; -7; 4; 5;
    % 1] etc.]. This would indicate that vertex 7 intrudes between vertices 3
    % and 4 in the resulting polygon. As another example, assume that vertex 2
    % and 7 both intersect simultaneously into other edges, between 9 and 10
    % for vertex 2, and between 4 and 5 for vertex 7. Then the
    % updatedVertexSequence will have values such as [1; 2; 3; -7; 4; 5; 6; 7;
    % 8; 9; -2; 10; 1]. Finish by flagging the intersection points as treated.
    %


    % STEP 3) Create a cell array to store which vertices will go into which
    % polytopes, one cell array for each polytope. Starting from the first
    % nonConvex vertex, proceed along upward in positive value until
    % returning back to start, jumping at any negative values to the
    % corresponding positive value, until reaching either positive or
    % negative valued version of the starting point. Then repeat this for
    % the negative valued version of the vertex ID.
    %
    % For example, consider the updatedVertexSequence as follows:
    %
    % originalSequence      = [1; 2; 3; 4; 5; 6; 7; 8; 9; 10; 11; 1]
    % flag_verticesToMerge  = [0; 0; 0; 0; 0; 0; 0; 0; 0;  1;  1; 0]  
    % updatedVertexSequence = [1; 2; 3; -7; 4; 5; 6; 7; 8; 9; -2; 10; 11; 1].
    % 
    % Vertex 2:
    % (forward) [2; 3; -7; 8; 9; -2] --> [2; 3; 7; 8; 9; 2]
    % (reverse) [-2; 10; 11; 1; 2];  --> [1; 2; 10; 11; 1]
    % Vertex 7:
    % (forward) [7; 8; 9; -2; 3; -7] --> [2; 3; 7; 8; 9; 2]
    % (reverse) [-7; 4; 5; 6; 7] --> [4; 5; 6; 7; 4]
    %
    % Another example:
    % originalSequence      = [1; 2; 3; 4; 5; 6; 7; 1]
    % flag_verticesToMerge  = [1; 1; 0; 0; 1; 0; 0; 1]
    % updatedVertexSequence = [1; -5; 2; 3; 4; 5; 6; 7; 1]
    % 
    % Vertex 5:
    %     (forward) [5; 6; 7; 1; -5]
    %     (reverse) [-5; 2; 3; 4; 5]

    % Step 4) Eliminate repeated points. Eliminate repeated polys.

    % Initialize output
    vertexSkeletonStartingPolytopes = struct;

    %%%%%%
    % STEP 1) find mergers
    mergedVertexIDList = fcn_INTERNAL_findRepeatedIntersectionsWithSameLocations(intersection_points, indices_repeated);

    %%%%%%
    % STEP 2) find updatedVertexSequence
    updatedVertexSequence = fcn_INTERNAL_updateVertexSequence(intersection_points, indices_repeated, flag_vertexIsNonConvex, boundaryEngagedAtMinCut);

    %%%%%
    % Step 3) 
    polytopeVertexIndices = fcn_INTERNAL_separateVertexSequenceIntoPolytopes(updatedVertexSequence);
       
    %%%%%
    % Step 4) 
    polytopeVertexIndicesNoRepeats = fcn_INTERNAL_removeRepeats(polytopeVertexIndices, mergedVertexIDList);

    % Output final results
    for ith_polytope = 1:length(polytopeVertexIndicesNoRepeats)
        vertexSkeletonStartingPolytopes.polytope(ith_polytope).vertices = intersection_points(polytopeVertexIndicesNoRepeats{ith_polytope},:);
    end

else
    warning('on','backtrace');
    warning('A vector was given that has dimension: %.0d, where 2D was expected. An error will be thrown here.',dimension_of_points);
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
        ([]), ...  % unit_normal_vectors
        ([]), ...  % unit_vertex_projection_vectors
        ([]), ... % vector_direction_of_unit_cut
        ([]),...  % flag_vertexIsNonConvex
        (0),...  % flag_plotEdgeGhostlines
        (0),...  % flag_plotVertexProjectionGhostlines
        ([]),... % plot_formatting
        (fig_num));  % fig_num

    colors = get(gca,'ColorOrder');

    for ith_polytope = 1:length(polytopeVertexIndicesNoRepeats)
        % Get the vertices for this polytope
        vertexPoints = vertexSkeletonStartingPolytopes.polytope(ith_polytope).vertices;

        % Set the color for this plot
        this_color_index = mod(ith_polytope,length(colors(:,1)))+1;
        plot_formatting.vertices_plot.Color = colors(this_color_index,:);
        plot_formatting.vertices_plot.vertexLabelsColor = colors(this_color_index,:);

        % Call the plotting function
        fcn_VSkel_plotPolytopeDetails(...
            vertexPoints,...
            ([]), ...  % unit_normal_vectors
            ([]), ...  % unit_vertex_projection_vectors
            ([]), ... % vector_direction_of_unit_cut
            ([]),...  % flag_vertexIsNonConvex
            (0),...  % flag_plotEdgeGhostlines
            (0),...  % flag_plotVertexProjectionGhostlines
            (plot_formatting),... % plot_formatting
            (fig_num));  % fig_num

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
function indices_in_this_merge = fcn_INTERNAL_findIndicesInThisIntersection(intersection_points, indices_repeated, this_index)
% Given a list of vertices, a list of possible repeats, and this index,
% checks to see if the repeated points are the SAME point, within a
% tolerance. For those that are the same, returns the indices that are all
% the same.

% Note the possible repeats
possibleRepeatVertices = intersection_points(indices_repeated,:);

% Check if any other intersection points are the same
% NOTE: may need to add a tolerance here
current_point = intersection_points(this_index,:);
differences_with_current = sum((possibleRepeatVertices - current_point).^2,2);
tolerance = 1E-6;
flag_vertices_same = (differences_with_current<tolerance);

% Find which vertices are intersecting here
indices_in_this_merge = indices_repeated(flag_vertices_same);

end % Ends fcn_INTERNAL_findIndicesInThisIntersection

%% fcn_INTERNAL_findRepeatedIntersectionsWithSameLocations
function mergedVertexIDList = fcn_INTERNAL_findRepeatedIntersectionsWithSameLocations(intersection_points, indices_repeated)
%%%%%%
% STEP 1) Loop through all indices_repeated, making sure not to analyze
% indicies previously treated, find all points in this intersection.
% Check if all these intersection points are all same. If so, flag
% that vertices are all merged into one. Save result as
%
%    flag_verticesToMerge
%
% with 0 meaning no mergers, 1 meaning the first merger, 2 meaning the
% second merger,
% etc. For example, flag_verticesToMerge would be, if vertex 4 is the
% first merger and is merging with 5 and 6 but not 7, 8, 9, etc. then
% flag_verticesToMerge would be [0 0 0 1 1 1 0 0 ...]. Finish by
% flagging the merged points as treated.

NumUniqueVerticies = length(intersection_points(:,1))-1;
% Initialize output
mergedVertexIDList = [(1:NumUniqueVerticies)'; 1];

% Initialize which vertices have been eliminated
indices_repeated_alreadyAnalyzed = indices_repeated*0;
flag_verticesToMerge = zeros(NumUniqueVerticies+1,1);
Num_mergers = 0; % Initialize merge count

for ith_repeatTest = 1:length(indices_repeated)
    % Make sure not to analyze indices that were previously analyzed
    if indices_repeated_alreadyAnalyzed(ith_repeatTest)==0
        this_index = indices_repeated(ith_repeatTest);
        indices_in_this_merge = fcn_INTERNAL_findIndicesInThisIntersection(intersection_points, indices_repeated, this_index);

        % Set indicies in for loop ahead of this point to not repeat
        % the search
        overlap =intersect(indices_in_this_merge,indices_repeated);
        indices_repeated_alreadyAnalyzed(find(overlap)) = 1; %#ok<FNDSB>

        % Increment number of mergers
        Num_mergers = Num_mergers+1;

        % Update the outputs
        flag_verticesToMerge(indices_in_this_merge,:) = Num_mergers;
        mergedVertexIDList(indices_in_this_merge) = min(indices_in_this_merge);
    end % Ends test to see if this is a repeated analysis
end % Ends for loop through repeated indices
end % Ends fcn_INTERNAL_findRepeatedIntersectionsWithSameLocations


%% fcn_INTERNAL_updateVertexSequence
function updatedVertexSequence = fcn_INTERNAL_updateVertexSequence(intersection_points, indices_repeated, flag_vertexIsNonConvex, boundaryEngagedAtMinCut)
%%%%%%%%%%%%%%%%%%%%
% STEP 2) Of indicies_repeated, making sure not to analyze non-convex
% indicies previously treated, loop through each one, finding insertion
% points for each indicating where polytopes are separating. New polytopes
% must be made at each of these vertices. This means that verticies have to
% be grouped into sub-polytopes in a later step. This is difficult, because
% the grouping of points into polytopes means all indicies - even ones not
% yet processed - need to be accounted for, otherwise polytopes can be
% created that would then have to themselves be broken again in later
% steps. Thus, the point of this step is to process the insertion points of
% all vertices, thus defining where polytope "breaks" occur from one
% polytope to another.  This is done by creating and updating a list of
% point indicies and insert a -N into the list at the location of break
% point, where N indicates the location where nonConvex vertex N intrudes.
% So the updatedVertexSequence will have values such as [1; 2; 3; -7; 4; 5;
% 1] etc.]. This would indicate that vertex 7 intrudes between vertices 3
% and 4 in the resulting polygon. As another example, assume that vertex 2
% and 7 both intersect simultaneously into other edges, between 9 and 10
% for vertex 2, and between 4 and 5 for vertex 7. Then the
% updatedVertexSequence will have values such as [1; 2; 3; -7; 4; 5; 6; 7;
% 8; 9; -2; 10; 1]. Finish by flagging the intersection points as treated.


NumUniqueVerticies = length(intersection_points(:,1))-1;

% Initialize output
updatedVertexSequence = [(1:NumUniqueVerticies)'; 1];

% Initialize which vertices have been eliminated
indices_repeated_alreadyAnalyzed = indices_repeated*0;

% Loop through each indices_repeated
for ith_intersection = 1:length(indices_repeated)
    % Make sure not to analyze indices that were previously analyzed
    if indices_repeated_alreadyAnalyzed(ith_intersection)==0
        this_index = indices_repeated(ith_intersection);
        indices_in_this_merge = fcn_INTERNAL_findIndicesInThisIntersection(intersection_points, indices_repeated, this_index);

        % Set indicies in for loop ahead of this point to not repeat
        % the search
        overlap =intersect(indices_in_this_merge,indices_repeated);
        indices_repeated_alreadyAnalyzed(find(overlap)) = 1; %#ok<FNDSB>

        % Confirm at least one of thse is NOT convex
        flag_thisMergeIsNOTAllConvex = any(flag_vertexIsNonConvex(indices_in_this_merge)==1);

        if flag_thisMergeIsNOTAllConvex
            % Find which boundary insertion is being done upon
            insertionPoint = boundaryEngagedAtMinCut(ith_intersection);
            
            % Find insertion point within updatedVertexSequence
            insertionStart = find(updatedVertexSequence==insertionPoint,1,'first');

            % Perform insertion
            updatedVertexSequence = [updatedVertexSequence(1:insertionStart,:); -1*this_index; updatedVertexSequence(insertionStart+1:end)];

        end
    end % Ends test to see if this is a repeated analysis
end % Ends for loop through repeated indices

% For debugging
if 1==0
    disp(updatedVertexSequence);
end
end % Ends fcn_INTERNAL_updateVertexSequence

%% fcn_INTERNAL_separateVertexSequenceIntoPolytopes
function polytopeVertexIndices = fcn_INTERNAL_separateVertexSequenceIntoPolytopes(updatedVertexSequence)
% STEP 3) Create a cell array to store which vertices will go into which
% polytopes, one cell array for each polytope. Starting from the first
% nonConvex vertex, proceed along upward in positive value until
% returning back to start, jumping at any negative values to the
% corresponding positive value, until reaching either positive or
% negative valued version of the starting point. Then repeat this for
% the negative valued version of the vertex ID.
%
% For example, consider the updatedVertexSequence as follows:
%
% originalSequence      = [1; 2; 3; 4; 5; 6; 7; 8; 9; 10; 11; 1]
% flag_verticesToMerge  = [0; 0; 0; 0; 0; 0; 0; 0; 0;  1;  1; 0]
% updatedVertexSequence = [1; 2; 3; -7; 4; 5; 6; 7; 8; 9; -2; 10; 11; 1].
%
% Vertex 2:
% (forward) [2; 3; -7; 8; 9; -2] --> [2; 3; 7; 8; 9; 2]
% (reverse) [-2; 10; 11; 1; 2];  --> [1; 2; 10; 11; 1]
% Vertex 7:
% (forward) [7; 8; 9; -2; 3; -7] --> [2; 3; 7; 8; 9; 2]
% (reverse) [-7; 4; 5; 6; 7] --> [4; 5; 6; 7; 4]
%
% Another example:
% originalSequence      = [1; 2; 3; 4; 5; 6; 7; 1]
% flag_verticesToMerge  = [1; 1; 0; 0; 1; 0; 0; 1]
% updatedVertexSequence = [1; -5; 2; 3; 4; 5; 6; 7; 1]
%
% Vertex 5:
%     (forward) [5; 6; 7; 1; -5]
%     (reverse) [-5; 2; 3; 4; 5]

negativeIndices = updatedVertexSequence(updatedVertexSequence<0);
numNegatives = length(negativeIndices);

numPolytopes = max(1,numNegatives*2);

% Initialize the output
polytopeVertexIndices = cell(numPolytopes,1);

if numNegatives==0
    polytopeVertexIndices{1} = updatedVertexSequence;
else
    for ith_polytope = 1:numNegatives
        this_negative_index = negativeIndices(ith_polytope);

        % Find the negative valued sequence
        negative_sequence = fcn_INTERNAL_extractSequence(updatedVertexSequence,this_negative_index);
        polytopeVertexIndices{ith_polytope*2 - 1} = negative_sequence;

        % Find the positive valued sequence
        positive_sequence = fcn_INTERNAL_extractSequence(updatedVertexSequence,this_negative_index*-1);
        polytopeVertexIndices{ith_polytope*2} = positive_sequence;
    end

end

end % Ends fcn_INTERNAL_separateVertexSequenceIntoPolytopes


%% fcn_INTERNAL_extractSequence
function output_sequence_sorted = fcn_INTERNAL_extractSequence(inputSequence,indexToFind)

% For debugging
if 1==0
    inputSequence = [1; 2; 3; -7; 4; 5; 6; 7; 8; 9; -2; 10; 11; 1];
    indexToFind = 7;
end

starting_point = fcn_INTERNAL_findLocationInSequence(inputSequence, indexToFind);
NinputSequence = length(inputSequence);

% Initialize variables
NumVerticesSearched = 1; % Incrementa a count of verticies, to avoid infinite while loops
flag_keepGoing = 1; % Flag that keeps us in while loop
temp_output_sequence(1,1) = indexToFind; % Initialize the output sequence
previousLocation = starting_point;

while(flag_keepGoing)
    NumVerticesSearched = NumVerticesSearched+1;

    nextLocation = previousLocation+1;

    % Check for wrap-around
    if nextLocation==NinputSequence
        nextLocation = 1;
    end
    nextIndex = inputSequence(nextLocation);
    temp_output_sequence(NumVerticesSearched,1) = nextIndex;

    if abs(nextIndex)==abs(indexToFind)
        flag_keepGoing = 0;
    end
    if nextIndex<0
        previousLocation = fcn_INTERNAL_findLocationInSequence(inputSequence, -1*nextIndex);
    else
        previousLocation = nextLocation;
    end

    if NumVerticesSearched>NinputSequence
        fprintf(1,'Input sequence: \t');
        disp(inputSequence');
        fprintf(1,'Current output sequence: \t');
        disp(temp_output_sequence');
        warning('on','backtrace');
        warning('When searching the above input sequence for the indexToFind: %.0f, ended adding more vertices than the original sequence. Code seems to be looping inescapably. An error will be thrown now.', indexToFind);
        error('Too many vertices encountered in constructing polytope');
    end

end % ends while loop

% Keep only the positive values
output_sequence = abs(temp_output_sequence);

% Sort the sequence
output_sequence_sorted = fcn_INTERNAL_sortVertexSequence(output_sequence);

end % Ends fcn_INTERNAL_extractSequence

%% fcn_INTERNAL_removeRepeats
function polytopeVertexIndicesNoRepeats = fcn_INTERNAL_removeRepeats(polytopeVertexIndices, mergedVertexIDList)
% Elminate merged points and repeated polytopes

% Used mergedVertexIDList as a dictionary to indicate repeated indices
polytopeVertexIndicesNoMergeRepeats = cell(length(polytopeVertexIndices),1);
longestVertexList = 0;
for ith_polytope = 1:length(polytopeVertexIndices)
    oldIndexList = polytopeVertexIndices{ith_polytope};
    newIndexList = mergedVertexIDList(oldIndexList); % apply the dictionary
    newIndexListUnique = unique(newIndexList,'stable'); % remove repeats
    newPolytopeIndexList = [newIndexListUnique; newIndexListUnique(1,1)]; % % must reclose the loop, since the unique function deletes repeats
    polytopeVertexIndicesNoMergeRepeats{ith_polytope} = newPolytopeIndexList;

    % Update length of longest vertex list
    longestVertexList = max(longestVertexList,length(newPolytopeIndexList));
end

% Make sure there are no polynomial repeats - this can happen when multiple
% vertices intersect at same time (symmetry)
verticesInMatrixForm = zeros(length(polytopeVertexIndices),longestVertexList);

for ith_polytope = 1:length(polytopeVertexIndices)
    indexList = polytopeVertexIndicesNoMergeRepeats{ith_polytope};
    verticesInMatrixForm(ith_polytope,1:length(indexList)) = indexList';    
end
uniqueRows = unique(verticesInMatrixForm,'rows','stable');

% Deal out only the unique rows
polytopeVertexIndicesNoRepeats = cell(length(uniqueRows(:,1)),1);
for ith_polytope = 1:length(uniqueRows(:,1))
    thisRow = uniqueRows(ith_polytope,:);
    polytopeVertexIndicesNoRepeats{ith_polytope} =thisRow(thisRow>0);
end


end % Ends fcn_INTERNAL_removeRepeats

%% fcn_INTERNAL_findLocationInSequence
function starting_point = fcn_INTERNAL_findLocationInSequence(inputSequence, indexToFind)
% Finds the index where the inputSequence has a value of indexToFind
starting_point = find(inputSequence(1:end-1,:)==indexToFind);
if isempty(starting_point) || length(starting_point)>1
    disp(inputSequence);
    warning('on','backtrace');
    warning('An input sequence was given, shown above, that contains %.0f of the indexToFind: %.0f. An error will be thrown now.',length(starting_point), indexToFind);
    error('Unexpected result: a vertex is either missing or was found twice in a vertex sequence.');
end
end % Ends fcn_INTERNAL_findLocationInSequence


%% fcn_INTERNAL_sortVertexSequence

function sortedSequence = fcn_INTERNAL_sortVertexSequence(inputSequence)
% Rearranges a vertex sequence so that the lowest vertex number always
% starts first

min_vertex_index = min(inputSequence);
starting_point = find(inputSequence(1:end-1,:)==min_vertex_index);
if isempty(starting_point) || length(starting_point)>1
    disp(inputSequence);
    warning('on','backtrace');
    warning('An input sequence was given, shown above, that seems to contain either 0 or more than 1 minimum vertex: %.0f.',starting_point);
    error('Unexpected result: a vertex is either missing or was found twice in a vertex sequence.');
end
sortedSequence = [inputSequence(starting_point:end-1,:); inputSequence(1:starting_point-1); inputSequence(starting_point,:)];
end % Ends fcn_INTERNAL_sortVertexSequence

