function [unit_normal_vectors, ...
    vertex_projection_vectors] = ...
    fcn_VSkel_fcn_polytopeFindUnitDirectionVectors(vertices,varargin)

%% fcn_VSkel_fcn_polytopeFindUnitDirectionVectors
% finds the vector_direction_of_unit_cut to use out of each vertex point,
% e.g. the direction and distance needed to move each point, given a
% unit edge cut
%
% FORMAT:
%
% [unit_normal_vectors, ...
%     vertex_projection_vectors] = ...
%     fcn_VSkel_fcn_polytopeFindUnitDirectionVectors(vertices, (fig_num))
%
% INPUTS:
%
%     vertices: a (M+1)-by-2 matrix of xy points with row1 = rowm+1, where
%         M is the number of the individual polytope vertices
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
%     unit_normal_vectors: a cell array of dimension M, where
%     each index 1:M stores a N x 2 array of the unit vectors that point
%     inward as measured from one vertex to the next.
%
%     vertex_projection_vectors: a cell array of M, where each index 1:M
%     stores a N x 2 array of the unit vectors that point
%     away from the vertices into the nested shape inside, with M = 1 being
%     the starting unit vectors and N being smaller and smaller for each M value.
%
% DEPENDENCIES:
%
%     fcn_DebugTools_checkInputsToFunctions
%
% % EXAMPLES:
% For additional examples, see: script_test_fcn_VSkel_fcn_polytopeFindUnitDirectionVectors
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

if flag_check_inputs
    % Are there the right number of inputs?
    if nargin < 1 || nargin > 2
        error('Incorrect number of input arguments')
    end

    % Check the vertices input
    fcn_DebugTools_checkInputsToFunctions(...
        vertices, '2or3column_of_numbers');

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

% Is this 2D or 3D
dimension_of_points = length(vertices(1,:));

% Calculate the unit vectors for each edge
if 2==dimension_of_points
    % find distances and unit vectors from vertex to vertex
    difference_vertex_to_vertex = vertices(2:end,:)-vertices(1:end-1,:);
    distances_vertex_to_vertex = sum(difference_vertex_to_vertex.^2,2).^0.5;
    unit_vectors_vertex_to_vertex = difference_vertex_to_vertex./distances_vertex_to_vertex;
    
    % Repeat last vector so it has same length as points
    unit_vectors_vertex_to_vertex = [unit_vectors_vertex_to_vertex; unit_vectors_vertex_to_vertex(1,:)];

    % Rotate by 90 degrees to get unit normal vectors
    unit_normal_vectors = unit_vectors_vertex_to_vertex*[0 1; -1 0];

    % Find the vertex vectors - do the calculation for vector 2 to end
    pseudo_vertex_projection_vectors = (unit_normal_vectors(2:end,:) + unit_normal_vectors(1:end-1,:))/2;
    pseudo_vertex_projection_vectors = [pseudo_vertex_projection_vectors(end,:);pseudo_vertex_projection_vectors];  % vector 1 is the same as the last one

    pseudo_vertex_projection_vector_lengths = sum(pseudo_vertex_projection_vectors.^2,2).^0.5;
    vertex_projection_vectors = pseudo_vertex_projection_vectors./pseudo_vertex_projection_vector_lengths;
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


    grid on
    grid minor
    hold on
    axis equal


    % Find size of vertex domain
    size = max(max(vertices)) - min(min(vertices));
    nudge = size*0.006;

    % Find the modpoints for each vertex
    midpoints = (vertices(2:end,:)+vertices(1:end-1,:))/2;
    midpoints = [midpoints; midpoints(1,:)]; % Repeat first row, to last, to match how point is similarly repeated

    % Plot the polytope in black dots connected by lines
    plot(vertices(:,1),vertices(:,2),'b.-','Linewidth',2, 'MarkerSize',10);

    % Label the vertices with their numbers
    for ith_vertex = 1:length(vertices(:,1))-1
        text(vertices(ith_vertex,1)+nudge,vertices(ith_vertex,2),...
            sprintf('%.0d',ith_vertex));
    end

    % Draw the unit vectors
    quiver(midpoints(1:end-1,1),midpoints(1:end-1,2),unit_normal_vectors(1:end-1,1),unit_normal_vectors(1:end-1,2),0,'r');

    % Draw the vertex_projection_vectors 
    quiver(vertices(1:end-1,1),vertices(1:end-1,2), vertex_projection_vectors(1:end-1,1),vertex_projection_vectors(1:end-1,2),0,'g');

    % Make axis slightly larger?
    if flag_rescale_axis
        temp = axis;
        %     temp = [min(points(:,1)) max(points(:,1)) min(points(:,2)) max(points(:,2))];
        axis_range_x = temp(2)-temp(1);
        axis_range_y = temp(4)-temp(3);
        percent_larger = 0.3;
        axis([temp(1)-percent_larger*axis_range_x, temp(2)+percent_larger*axis_range_x,  temp(3)-percent_larger*axis_range_y, temp(4)+percent_larger*axis_range_y]);
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