function [cut_distance, vertexSkeleton] = ...
    fcn_VSkel_polytopeFindVertexSkeleton(vertices, varargin)
% Calculates the VertexSkeleton for a polytope, i.e. where the vertices
% would land if the polytope were shrunk.
%
% FORMAT:
%
% [cut_distance, vertexSkeleton] = ...
% fcn_VSkel_polytopeFindVertexSkeleton(vertices, varargin)
%
% INPUTS:
%
%     vertices: a m+1-by-2 matrix of xy points with row1 = rowm+1, where m is
%         the number of the individual polytope vertices
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
%     cut_distance: an array of M x 1, starting from 0 for M(1) to the
%     maximum cut distance that can be used, at M(end)
%
%     vertexSkeleton: a listing of the vertices for each cut distance,
%     using the following structure:
%
%        vertexSkeleton(depth_index).FIELDS
%
%     where depth_index refers to the index of the depth of cut. For
%     example, if the cut_distance is [0; 0.4; 0.7; 1], and the vertex
%     skeleton is requested for a cut distance of 0.2, then 0.2 lies
%     between the first depth (0) and the second depth (0.4). It will
%     correspond to the smallest of these (0), which is at index 1.
%     Therefore, the depth_index for these vertices is 1.
%
%     The subfields are indexed by the polytope_index, which counts the
%     number of active polytopes at this given cut depth:
%
%        vertexSkeleton(depth_index).polytope(polytope_index).SUBFIELDS
% 
%     where required subfields are as follows:
%
%        vertexSkeleton(depth_index).polytope(polytope_index).vertices                       = [1 2; 3 4; 5 6; 1 2]; % for 2D, the 1st and last points are same, so this is of size N+1 if there are N unique vertices 
%        vertexSkeleton(depth_index).polytope(polytope_index).vector_direction_of_unit_cut   = [1 2; 3 4; 5 6]       % same number of rows as vertices;
%
%     the following are also generated, but are optional:
%        vertexSkeleton(depth_index).polytope(polytope_index).unit_vertex_projection_vectors = [1 2; 3 4; 5 6];      % same number of rows as vertices; 
%        vertexSkeleton(depth_index).polytope(polytope_index).flag_vertexIsNonConvex         = [0; 0; 0; 0];         % same number of rows as vertices;
%        vertexSkeleton(depth_index).polytope(polytope_index).boundaries{boundary_index}     = [1 2];                % boundaries are defined, row-wise, as sequences of vertex IDs (rows in vertices list). Here, vertex 1 is connected to 2 as boundary 1. In 3D, the number of rows can be large but must start/end on same last point .
%        vertexSkeleton(depth_index).polytope(polytope_index).unit_normal_vector{edge_index} = [1 2];                % normal vector to the boundary, so will be a [1xD] where D is the dimension;
%
% TO DO: UPDATE:
% vertexSkeleton(depth_index).polytope(polytope_index).unit_vertex_projection_vectors = unit_vertex_projection_vectors;
% vertexSkeleton(depth_index).polytope(polytope_index).flag_vertexIsNonConvex         = flag_vertexIsNonConvex;
% vertexSkeleton(depth_index).polytope(polytope_index).boundaries{boundary_index}     = [];  % which vertices form each boundary
% vertexSkeleton(depth_index).polytope(polytope_index).boundaryEngagedAtMinCut        = boundaryEngagedAtMinCut;   % For each boundary, the unit normal vectors
% vertexSkeleton(depth_index).polytope(polytope_index).indices_repeated               = indices_repeated;   % For each boundary, the unit normal vectors
% vertexSkeleton(depth_index).polytope(polytope_index).intersection_points            = intersection_points;   % For each boundary, the unit normal vectors
%
%
% DEPENDENCIES:
%
%     fcn_DebugTools_checkInputsToFunctions
%     INTERNAL_fcn_findUnitDirectionVectors
%     fcn_VSkel_polytopeFindVertexAngles
%     fcn_VSkel_plotVertexSkeleton
%     
%
% % EXAMPLES:
% For additional examples, see: script_test_fcn_VSkel_polytopeFindVertexSkeleton
%
% This function was written on 2022_02_15 by S. Brennan
% Questions or comments? sbrennan@psu.edu
%

% Revision History:
% 2022_02_13 - S. Brennan
% -- first write of code
% -- pulled the function out of edge shrinking code
% 2025_04_25 by Sean Brennan
% -- added global debugging options
% -- switched input checking to fcn_DebugTools_checkInputsToFunctions


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

flag_do_debug = 1;

if flag_do_debug
    st = dbstack; %#ok<*UNRCH>
    fprintf(1,'STARTING function: %s, in file: %s\n',st(1).name,st(1).file);
    debug_fig_num = 999978; 
else
    debug_fig_num = -1; 
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
if (0==flag_max_speed)
    if flag_check_inputs
        % Are there the right number of inputs?
        if nargin < 1 || nargin > 2
            error('Incorrect number of input arguments')
        end

        % Check the vertices input
        % fcn_DebugTools_checkInputsToFunctions(...
        %     vertices, '2column_of_numbers');

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

% Did the user pass one polytope, or many? If just one, convert to cell
% array
if iscell(vertices)
    working_polytopes = vertices;
else
    working_polytopes = {vertices};  % The vertices we are now using are the starting vertices
end

% Initialize results
depth_index = 1; % Initialize iterations to 1
flag_stop_loop = 0;  % Set stop flag to 0 (e.g. "keep going!")
total_cut = 0;  % The total cut distance is set to zero
vertexSkeleton = struct;
vertexSkeleton(1).polytope(1) = struct;
vertexSkeleton(1).polytope(1).vertices = []; 


% Loop until the stop flag is set, e.g. loop "inward" to calculate the
% skeleton
while 0 == flag_stop_loop
    Npolytopes = length(working_polytopes);

    % Set flags to zero that are used to indicate if the polytope has been
    % processed
    flagsPolytopeIsDone = zeros(Npolytopes,1);
    maxCutDistancesEachPolytope = inf(Npolytopes,1);

    % Save the current iteration's vertices and cut_distance
    for polytope_index = 1:Npolytopes
        vertexSkeleton(depth_index).polytope(polytope_index).vertices{polytope_index} = working_polytopes{polytope_index};
    end
    cut_distance(depth_index,1) = total_cut; %#ok<AGROW>


    % Loop through the polytopes, getting projection vectors and finding
    % maximum cut depth

    for polytope_index = 1:Npolytopes
        thisPolytope = working_polytopes{polytope_index};

        %%%%%%
        % Find the unit vectors that point inward from each vertex point.
        vertices = thisPolytope;
        error('replace the function below');
        % [unit_normal_vectors, unit_vertex_projection_vectors, vector_direction_of_unit_cut, flag_vertexIsNonConvex]  = fcn_VSkel_polytopeFindUnitDirectionVectors(vertices,debug_fig_num);

        % Save the results into the projection vectors for this cut
        vertexSkeleton(depth_index).polytope(polytope_index).vertices = vertices;
        vertexSkeleton(depth_index).polytope(polytope_index).vector_direction_of_unit_cut   = vector_direction_of_unit_cut;
        

        % Check how many vertices we have. If only 2, then this polytope is
        % done shrinking
        Nvertices = length(thisPolytope(:,1));

        if 1>=Nvertices
            warning('on','backtrace');
            warning('A polynomial was given that does not have enough verticies. The row count was: %.0f, and 2 or more rows is required. An error will be thrown here.',Nvertices);
            error('Cannot shrink polytopes with ill-defined verticies.');
        elseif 2==Nvertices  % If 2, then no way to project, since this is the same point - remove this polytope from future calculations by NOT including into polytopesRemaining
            flagsPolytopeIsDone(polytope_index,1) = 1;
            maxCutDistancesEachPolytope(polytope_index,1) = 0;
            flag_vertexIsNonConvex = zeros(2,1);
            boundaryEngagedAtMinCut = zeros(2,1);
            indices_repeated = zeros(2,1);
            intersection_points = zeros(2,2);
            min_cut = 0;
        elseif 3==Nvertices % This is a simple line segment with just 2 points
            % Calculate the final cut distance using the first row to represent
            % the last point.
            last_cut_dist = sum((mean(vertices(1:2,:)) - vertices(1,:)).^2,2).^0.5; % Calculate this last cut distance
            maxCutDistancesEachPolytope(polytope_index,1) = last_cut_dist;
            flag_vertexIsNonConvex = zeros(3,1);
            boundaryEngagedAtMinCut = zeros(3,1);
            indices_repeated = zeros(3,1);
            intersection_points = zeros(3,2);
            min_cut = last_cut_dist;
        
        else % More then 3 vertices
            % Calculate intersection points
            max_edge_cuts = fcn_VSkel_polytopeFindMaxEdgeCut(vertices, unit_normal_vectors, unit_vertex_projection_vectors, (debug_fig_num));
            [sphereRadii, definingBoundaries] = fcn_VSkel_polytopeFindEnclosedSpheres(vertices, unit_normal_vectors, unit_vertex_projection_vectors, vector_direction_of_unit_cut, flag_vertexIsNonConvex, max_edge_cuts, (debug_fig_num));
            [min_cut, boundaryEngagedAtMinCut, indices_repeated, intersection_points] = ...
                fcn_VSkel_polytopeFindMinimumEnclosedSphere(vertices, sphereRadii, definingBoundaries, unit_normal_vectors, unit_vertex_projection_vectors, vector_direction_of_unit_cut, (debug_fig_num)); 
            maxCutDistancesEachPolytope(polytope_index,1) = min_cut;
        end

        % Save results for later use
        vertexSkeleton(depth_index).polytope(polytope_index).unit_vertex_projection_vectors = unit_vertex_projection_vectors;
        vertexSkeleton(depth_index).polytope(polytope_index).flag_vertexIsNonConvex         = flag_vertexIsNonConvex;
        vertexSkeleton(depth_index).polytope(polytope_index).boundaries{1}                  = [(1:length(vertices(:,1)))'; 1];  % which vertices form each boundary
        vertexSkeleton(depth_index).polytope(polytope_index).boundaryEngagedAtMinCut        = boundaryEngagedAtMinCut;   % For each boundary, the unit normal vectors
        vertexSkeleton(depth_index).polytope(polytope_index).indices_repeated               = indices_repeated;   % For each boundary, the unit normal vectors
        vertexSkeleton(depth_index).polytope(polytope_index).intersection_points            = intersection_points;   % For each boundary, the unit normal vectors
        vertexSkeleton(depth_index).polytope(polytope_index).min_cut                        = min_cut;   % For each boundary, the unit normal vectors
    end
    
    % Did all the polytopes disappear?
    if all(flagsPolytopeIsDone==1)
        flag_stop_loop = 1;
    else
        % No, they did not disappear, need to prep to shrink ones that
        % remain

        % Find which polytope is the limiting case. These are the "collapsed"
        % polytopes, e.g. the ones where vertices are re-arranging
        overall_minCut = min(maxCutDistancesEachPolytope(flagsPolytopeIsDone==0));

        indicesOfPolytopesCollapsing = abs(maxCutDistancesEachPolytope - overall_minCut)<1E-8;
        flagsPolytopeIsDone(indicesOfPolytopesCollapsing) = -1; % Flag the indicies of collapsing ones

        % Initialize polytopes remaining
        polytopesRemaining = cell(1,1);
        NpolytpoesRemaining = 0;

        % Loop through all the polytopes, moving vertices
        for polytope_index = 1:Npolytopes

            % Only need to move vertices if the polytope is NOT done. Polytopes
            % that are done are flagged as "1" values. "0" values mean
            % polytopes that have their vertices move, but none of the vertices
            % collapse. "-1" means the vertices collapse.
            if flagsPolytopeIsDone(polytope_index)<1
                vertices                     = vertexSkeleton(depth_index).polytope(polytope_index).vertices;
                vector_direction_of_unit_cut = vertexSkeleton(depth_index).polytope(polytope_index).vector_direction_of_unit_cut;
                flag_vertexIsNonConvex       = vertexSkeleton(depth_index).polytope(polytope_index).flag_vertexIsNonConvex;
                boundaryEngagedAtMinCut      = vertexSkeleton(depth_index).polytope(polytope_index).boundaryEngagedAtMinCut;
                indices_repeated             = vertexSkeleton(depth_index).polytope(polytope_index).indices_repeated ;
                intersection_points          = vertexSkeleton(depth_index).polytope(polytope_index).intersection_points;

                if flagsPolytopeIsDone(polytope_index)<0
                    % Find where the polytope vertices are moving for the "collapsed"
                    % polytopes
                    Nvertices = length(vertices(:,1));

                    if 3>Nvertices
                        warning('on','backtrace');
                        warning('A polynomial was given that does not have enough verticies. The row count was: %.0f, and 2 or more rows is required. An error will be thrown here.',Nvertices);
                        error('Cannot shrink polytopes with ill-defined verticies.');
                    elseif 3==Nvertices
                        vertexSkeletonStartingPolytopes.polytope(1).vertices = [mean(vertices(1:2,:),1); mean(vertices(1:2,:),1)];
                    else
                        % Note: the merger operation can take one polytope
                        % and break it into many polytopes. Hence, we need
                        % a for loop after this to bring the results back
                        % into the outer polytope list.
                        vertexSkeletonStartingPolytopes = ...
                            fcn_VSkel_polytopeMergeIntersectingVertices( ...
                            vertices, ...
                            flag_vertexIsNonConvex,...
                            boundaryEngagedAtMinCut, ...
                            indices_repeated, ...
                            intersection_points, ...
                            (debug_fig_num));
                    end

                    for ith_addedPolytope = 1:length(vertexSkeletonStartingPolytopes.polytope)
                        NpolytpoesRemaining = NpolytpoesRemaining + 1;
                        polytopesRemaining{NpolytpoesRemaining} = vertexSkeletonStartingPolytopes.polytope(ith_addedPolytope).vertices;
                    end
                else
                    % Need to move the vertices inward by the cut distance
                    NpolytpoesRemaining = NpolytpoesRemaining + 1;
                    polytopesRemaining{NpolytpoesRemaining} = vertices +  vector_direction_of_unit_cut* overall_minCut;

                end
            end

        end

        % Increment the depth_index count
        depth_index = depth_index+1;
        working_polytopes = polytopesRemaining;
        total_cut = total_cut + overall_minCut;
    end % Ends check to see if all flags are zero

end % Ends while loop

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
    fcn_VSkel_plotVertexSkeleton(cut_distance, vertexSkeleton, fig_num);
end % Ends flag_do_plot if statement

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file);
end

end % ends fucntion fcn_VSkel_polytopeFindVertexSkeleton


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

