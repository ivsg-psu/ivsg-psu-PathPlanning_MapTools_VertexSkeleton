function [unit_normal_vectors, unit_vertex_projection_vectors, vector_direction_of_unit_cut, flag_vertexIsNonConvex]  = ...
    fcn_VSkel_polytopeFindUnitDirectionVectors(vertices,varargin)

%% fcn_VSkel_polytopeFindUnitDirectionVectors
% finds the vector_direction_of_unit_cut to use out of each vertex point,
% e.g. the direction and distance needed to move each point, given a
% unit edge cut
%
% FORMAT:
%
% [unit_normal_vectors, ...
%     vector_direction_of_unit_cut, flag_vertexIsNonConvex] = ...
%     fcn_VSkel_polytopeFindUnitDirectionVectors(vertices, (fig_num))
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
%     unit_normal_vectors: an (M+1)-by-2 matrix of the unit vectors that
%     point inward as measured from one vertex to the next. The vector is
%     assumed to be attached to the start of the edge given by vertex M.
%
%     unit_vertex_projection_vectors: an (M+1)-by-2 matrix of the unit
%     vectors that project in the direction that each of the M verticies
%     will move. The vector is assumed to be attached to the start of the
%     edge given by vertex M.
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

        % Check the vertices input
        fcn_DebugTools_checkInputsToFunctions(...
            vertices, '2or3column_of_numbers');

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

% Is this 2D or 3D
dimension_of_points = length(vertices(1,:));
NumUniqueVerticies = length(vertices(:,1))-1;

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

    % Find the unit_vertex_projection_vectors
    if NumUniqueVerticies>2
        pseudo_vertex_projection_vectors = (unit_normal_vectors(2:end,:) + unit_normal_vectors(1:end-1,:))/2;
        pseudo_vertex_projection_vectors = [pseudo_vertex_projection_vectors(end,:);pseudo_vertex_projection_vectors];  % vector 1 is the same as the last one
        pseudo_vertex_projection_vector_lengths = sum(pseudo_vertex_projection_vectors.^2,2).^0.5;
        unit_vertex_projection_vectors = pseudo_vertex_projection_vectors./pseudo_vertex_projection_vector_lengths;
    elseif NumUniqueVerticies==2 % Line segment
        pseudo_vertex_projection_vectors = unit_vectors_vertex_to_vertex;
        pseudo_vertex_projection_vector_lengths = sum(pseudo_vertex_projection_vectors.^2,2).^0.5;
        unit_vertex_projection_vectors = pseudo_vertex_projection_vectors./pseudo_vertex_projection_vector_lengths;
    elseif NumUniqueVerticies==1 % Point
        pseudo_vertex_projection_vectors = zeros(size(unit_vectors_vertex_to_vertex));
        unit_vertex_projection_vectors = pseudo_vertex_projection_vectors;
    else
        warning('on','backtrace');
        warning('Expecting 1 or more points, but no vertices are available for vertex projection calculations! Throwing an error');
        error('Do not know how to calculate vertex projection for one point!')
    end



    % Find the vector direction of unit cuts
    % See the documentation.
    if NumUniqueVerticies>2
        vector_sums = sum(unit_normal_vectors.*unit_vertex_projection_vectors,2);
        d = 1./vector_sums;
    elseif NumUniqueVerticies==2 || NumUniqueVerticies==1 % Line segment or point
        d = 1;
    else
        warning('on','backtrace');
        warning('Expecting 2 or more points, but less than 2 are available for vertex projection calculations! Throwing an error');
        error('Do not know how to calculate vertex projection for one point!')
    end
    vector_direction_of_unit_cut = d.*unit_vertex_projection_vectors;

    % Check which ones are NOT convex. This is done by doing cross product of
    % vectors in sequence to see if their unit normals are in same or
    % opposite directions
    cross_products = cross([unit_normal_vectors(1:end-1,:) zeros(NumUniqueVerticies,1)],[unit_normal_vectors(2:end,:) zeros(NumUniqueVerticies,1)],2);
    cross_products = [cross_products(end,:); cross_products];
    cross_product_results = cross_products(:,3);
    flag_vertexIsNonConvex = cross_product_results<0;
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
    fcn_VSkel_plotPolytopeDetails(...
        vertices,...
        (unit_normal_vectors), ...
        (unit_vertex_projection_vectors), ...
        (vector_direction_of_unit_cut), ...
        (flag_vertexIsNonConvex),...
        (1),...  % flag_plotEdgeGhostlines
        (1),...  % flag_plotVertexProjectionGhostlines
        ([]),...  % plot_formatting
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