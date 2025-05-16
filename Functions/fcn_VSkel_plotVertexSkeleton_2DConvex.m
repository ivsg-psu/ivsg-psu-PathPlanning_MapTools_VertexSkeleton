function [h_fig] = ...
    fcn_VSkel_plotVertexSkeleton_2DConvex(vertices, projection_vectors, cut_distance, varargin)
% plots the skelton of a polytope, e.g. where the polytope will shrink if
% the edges are all brought in at exactly the same rate.
%
% FORMAT:
%
% fcn_VSkel_plotVertexSkeleton_2DConvex(vertices, projection_vectors, cut_distance, (fig_num))
%
% INPUTS:
%
%     vertices: a cell array of dimension M, where each index 1:M
%     stores a N x 2 array of the coordinates of the nested shape inside,
%     with M = 1 being the starting shape (and with dimension K+1, where K
%     is number of vertices given), and N being smaller and smaller for
%     each M value.
%
%     projection_vectors: a cell array of M, where each index 1:M
%     stores a N x 2 array of the unit vectors that point
%     away from the vertices of the nested shape inside, with M = 1 being
%     the starting unit vectors and N being smaller and smaller for each M value.
%
%     cut_distance: an array of 1 x M, starting from 0 for M(1) to the
%     maximum cut distance that can be used, at M(end)
%
%    (OPTIONAL INPUTS)
%
%     fig_num: a figure number to plot results.
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
% For additional examples, see: script_test_fcn_VSkel_plotVertexSkeleton_2DConvex
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
if (0==flag_max_speed)
    if flag_check_inputs
        % Are there the right number of inputs?
        narginchk(3,4);

        % Check the cut_distance input
        fcn_DebugTools_checkInputsToFunctions(...
            cut_distance, '1column_of_numbers');

    end
end


% Does user want to show the plots?
flag_do_plot = 0; % Default is no plotting
if  (4 == nargin) && (0==flag_max_speed) % Only create a figure if NOT maximizing speed
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

h_fig = figure(fig_num);
grid on
grid minor
hold on
axis equal

% Grab the original vertices
original_vertices = vertices{1};

% Plot the original polytope in red using the vertices
h_plot = plot(original_vertices(:,1),original_vertices(:,2),'r-','Linewidth',2);
last_color = get(h_plot,'Color');

% Find size of vertices so we can figure out nudging
size = max(max(original_vertices)) - min(min(original_vertices));
nudge = size*0.003;

% Number the original vertices with labels, nudging a litte so they don't
% land right on top of the points
for ith_vertex = 1:(length(original_vertices(:,1))-1)
    text(original_vertices(ith_vertex,1)+nudge,original_vertices(ith_vertex,2),...
        sprintf('%.0d',ith_vertex));
end

% Plot each contraction
Ncontractions = length(cut_distance);
for ith_contraction = 2:Ncontractions

    % Do calculations to determine vector start points, length of vectors,
    % and vectors themselves
    cut_length = cut_distance(ith_contraction)-cut_distance(ith_contraction-1);
    starting_points = vertices{ith_contraction-1}(1:end-1,:);
    vectors_from_starting_points = ...
        projection_vectors{ith_contraction-1}(1:end-1,:)*cut_length;

    % Plot the vectors as arrows going out from the starting point to
    % ending point. Use the color from the previous plot of the points, so
    % the vectors are same color as the points they originate from.
    quiver(starting_points(:,1),starting_points(:,2),vectors_from_starting_points(:,1),vectors_from_starting_points(:,2),0,'Color',last_color);

    % Plot the ending points, and save the color since they will be the new
    % start points on the next iteration.
    ending_points = vertices{ith_contraction}(1:end-1,:);
    h_plot = plot(ending_points(:,1),ending_points(:,2),'.','Markersize',20);
    last_color = get(h_plot,'Color');

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
    % Nothing to do here - it's a plotting function
end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file);
end

end % Ends fcn_VSkel_plotVertexSkeleton_2DConvex

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