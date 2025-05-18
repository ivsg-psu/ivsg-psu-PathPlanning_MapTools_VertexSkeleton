function [h_fig] = ...
    fcn_VSkel_plotVertexSkeleton(cut_distance, vertexSkeleton, varargin)
% plots the vertex skelton of a polytope, e.g. where the polytope will move if
% the edges are all cut with exactly the same rates.
%
% FORMAT:
%
% fcn_VSkel_plotVertexSkeleton(cut_distance, vertexSkeleton, (fig_num))
%
% INPUTS:
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
%     h_fig: a handle to the resulting figure
%
% DEPENDENCIES:
%
%     fcn_DebugTools_checkInputsToFunctions
%
% % EXAMPLES:
% For additional examples, see: script_test_fcn_VSkel_plotVertexSkeleton
%
% This function was written on 2022_02_15 by S. Brennan
% Questions or comments? sbrennan@psu.edu
%

% Revision History:
% 2022_02_16 - S. Brennan
% -- first write of code
% 2025_04_25 by Sean Brennan
% -- added global debugging options
% -- switched input checking to fcn_DebugTools_checkInputsToFunctions
% 2025_04_29 by Sean Brennan
% -- first written by S. Brennan by pulling function out of MapGen_polytopeFindVertexSkeleton
% 2025_05_15 by Sean Brennan
% -- modified to use the new vertexSkeleton structure

% TO DO
% -- none

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==3 && isequal(varargin{end},-1))
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
if (0==flag_max_speed)
    if flag_check_inputs
        % Are there the right number of inputs?
        narginchk(2,3);

        % Check the cut_distance input
        fcn_DebugTools_checkInputsToFunctions(...
            cut_distance, '1column_of_numbers');

    end
end


% Does user want to show the plots?
flag_do_plot = 1; % Default is for plotting
if  (3 == nargin) && (0==flag_max_speed) % Only create a figure if NOT maximizing speed
    temp = varargin{end}; % Last argument is always figure number
    if ~isempty(temp) % Make sure the user is not giving empty input
        fig_num = temp;
        flag_do_plot = 1; % Set flag to do plotting
    end
else
    fig = gcf; % Get current figure number
    fig_num = fig.Number; 

    if flag_do_debug % If in debug mode, do plotting but to an arbitrary figure number
        fig_for_debug = 454654; %#ok<NASGU>
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

% check whether the figure already has data
h_fig = figure(fig_num);
flag_rescale_axis = 0;
if isempty(get(h_fig,'Children'))
    flag_rescale_axis = 1;
end

grid on
grid minor
hold on
% axis equal



NoriginalPoloytopes = length(vertexSkeleton(1).polytope);

all_original_vertices = [];
for ith_polytope = 1:NoriginalPoloytopes
    % Grab the original vertices
    original_vertices = vertexSkeleton(1).polytope(ith_polytope).vertices;

    % Save vertices
    all_original_vertices = [all_original_vertices; original_vertices]; %#ok<AGROW>
end

% Find size of vertex domain
max_XY = max(all_original_vertices);
min_XY = min(all_original_vertices);
sizePlot = max(max_XY) - min(min_XY);
nudge = sizePlot*0.006;

% Make axis bigger
if flag_rescale_axis
    percent_larger = 0.3;

    axis_range_x = max_XY(1,1)-min_XY(1,1);
    axis_range_y = max_XY(1,2)-min_XY(1,2);

    if (0==axis_range_x)
        axis_range_x = 2/percent_larger;
    end
    if (0==axis_range_y)
        axis_range_y = 2/percent_larger;
    end
    axis([min_XY(1,1)-percent_larger*axis_range_x, max_XY(1,1)+percent_larger*axis_range_x,  min_XY(1,2)-percent_larger*axis_range_y, max_XY(1,2)+percent_larger*axis_range_y]);
end
goodAxis = axis;



% Plot each contraction
Ncontractions = length(cut_distance(:,1));
for ith_contraction = 1:Ncontractions

    % Do calculations to determine vector start points, length of vectors,
    % and vectors themselves
    if ith_contraction~=Ncontractions
        cut_length = cut_distance(ith_contraction+1)-cut_distance(ith_contraction);
    else
        cut_length = 0;
    end
    
    Npolytopes = length(vertexSkeleton(ith_contraction).polytope);
    last_color = zeros(Npolytopes,3);
    for ith_polytope = 1:Npolytopes

        % Grab the starting_points
        starting_points = vertexSkeleton(ith_contraction).polytope(ith_polytope).vertices;

        % Plot the starting_points
        h_plot = plot(starting_points(:,1),starting_points(:,2),'.-','Linewidth',0.5, 'MarkerSize',20);
        last_color(ith_polytope,:) = get(h_plot,'Color');


        % Number the original vertices with labels, nudging a litte so they don't
        % land right on top of the points
        for ith_vertex = 1:(length(starting_points(:,1))-1)
            text(starting_points(ith_vertex,1)+nudge,starting_points(ith_vertex,2),...
                sprintf('%.0d',ith_vertex), 'Color',last_color(ith_polytope,:));
        end


        vectors_from_starting_points = ...
            vertexSkeleton(ith_contraction).polytope(ith_polytope).vector_direction_of_unit_cut * cut_length;

        % Plot the vectors as arrows going out from the starting point to
        % ending point. Use the color from the previous plot of the points, so
        % the vectors are same color as the points they originate from.
        quiver(starting_points(:,1),starting_points(:,2),vectors_from_starting_points(:,1),vectors_from_starting_points(:,2),0,'Color',last_color(ith_polytope,:), 'LineWidth',3,'MaxHeadSize',0.05);
    end

end
axis(goodAxis);
axis equal

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
    % Nothing to do here - it's a plotting function
end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file);
end

end % Ends fcn_VSkel_plotVertexSkeleton

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