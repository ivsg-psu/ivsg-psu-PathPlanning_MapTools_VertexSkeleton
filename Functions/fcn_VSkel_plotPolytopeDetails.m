function h_fig = fcn_VSkel_plotPolytopeDetails(vertices, varargin)

%% fcn_VSkel_plotPolytopeDetails
% plots a single poloytpe allowing added user-defined details
%
% FORMAT:
%
% h_fig =  fcn_VSkel_plotPolytopeDetails(...
%        vertices,...
%        (unit_normal_vectors), ...  % unit_normal_vectors
%        (unit_vertex_projection_vectors), ...  % unit_vertex_projection_vectors
%        (vector_direction_of_unit_cut), ... % vector_direction_of_unit_cut
%        (flag_vertexIsNonConvex),...  % flag_vertexIsNonConvex
%        (flag_plotEdgeGhostlines),...  % flag_plotEdgeGhostlines
%        (flag_plotVertexProjectionGhostlines),...  % flag_plotVertexProjectionGhostlines
%        (plot_formatting),... % plot_formatting
%        (fig_num));  % fig_num
%
% INPUTS:
%
%     vertices: a (M+1)-by-2 matrix of xy points with row1 = rowm+1, where
%         M is the number of the individual polytope vertices
%
%     (OPTIONAL INPUTS)
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
%     flag_plotEdgeGhostlines: plots edge ghostlines showing where tangents
%     can occur
%
%     plot_formatting: a structure specifying the plot style. For any
%     fields left empty, defaults are used. The defaults are
%
%            plot_formatting.vertices_plot.style = 'b.-';
%            plot_formatting.vertices_plot.LineWidth = 2;
%            plot_formatting.vertices_plot.MarkerSize = 20;
%            plot_formatting.vertices_plot.Color = [0 0 1];
%            plot_formatting.vertices_plot.vertexLabelsColor = [0 1 0];
%            plot_formatting.vertices_plot.edgeLabelsColor = [0 0 1];
%            
%            plot_formatting.edgeGhostLines_plot.style = '-';
%            plot_formatting.edgeGhostLines_plot.LineWidth = 0.5;
%            plot_formatting.edgeGhostLines_plot.MarkerSize = 0.1;
%            plot_formatting.edgeGhostLines_plot.Color = 0.7*[1 1 1];
%            
%            plot_formatting.vertexProjectionGhostLines_plot.style = '-';
%            plot_formatting.vertexProjectionGhostLines_plot.LineWidth = 0.5;
%            plot_formatting.vertexProjectionGhostLines_plot.MarkerSize = 0.1;
%            plot_formatting.vertexProjectionGhostLines_plot.Color = 0.7*[0 1 0];
%            
%            plot_formatting.unitNormalVectors_plot.style = 'r';
%            plot_formatting.unitNormalVectors_plot.LineWidth = 0.5;
%            plot_formatting.unitNormalVectors_plot.MarkerSize = 0.1;
%            plot_formatting.unitNormalVectors_plot.Color = [1 0 0];
%            
%            plot_formatting.unitVertexProjectionVectors_plot.style = 'g';
%            plot_formatting.unitVertexProjectionVectors_plot.LineWidth = 3;
%            plot_formatting.unitVertexProjectionVectors_plot.MarkerSize = 0.1;
%            plot_formatting.unitVertexProjectionVectors_plot.Color = [0 1 0];
%            
%            plot_formatting.vectorDirectionOfUnitCut_plot.style = '-';
%            plot_formatting.vectorDirectionOfUnitCut_plot.LineWidth = 2;
%            plot_formatting.vectorDirectionOfUnitCut_plot.MarkerSize = 0.1;
%            plot_formatting.vectorDirectionOfUnitCut_plot.Color = [0 0.5 0];
%            
%            plot_formatting.vertexIsNonConvex_plot.style = '.';
%            plot_formatting.vertexIsNonConvex_plot.LineWidth = 2;
%            plot_formatting.vertexIsNonConvex_plot.MarkerSize = 20;
%            plot_formatting.vertexIsNonConvex_plot.Color = [1 0 0];
%
%     fig_num: a figure number to plot results. If set to -1, skips any
%     input checking or debugging, no figures will be generated, and sets
%     up code to maximize speed. As well, if given, this forces the
%     variable types to be displayed as output and as well makes the input
%     check process verbose.
%
% OUTPUTS:
%
%     h_fig: the figure handle
%
% DEPENDENCIES:
%
%     fcn_DebugTools_checkInputsToFunctions
%
% % EXAMPLES:
% For additional examples, see: script_test_fcn_VSkel_plotPolytopeDetails
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
if (nargin==9 && isequal(varargin{end},-1))
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
        narginchk(1,9);

        % Check the vertices input
        fcn_DebugTools_checkInputsToFunctions(...
            vertices, '2or3column_of_numbers');

    end
end

% Does user want to specify unit_normal_vectors?
unit_normal_vectors = [];
if nargin>=2
    temp = varargin{1};
    if ~isempty(temp)
        unit_normal_vectors = temp;
    end
end


% Does user want to specify unit_vertex_projection_vectors?
unit_vertex_projection_vectors = [];
if nargin>=3
    temp = varargin{2};
    if ~isempty(temp)
        unit_vertex_projection_vectors = temp;
    end
end


% Does user want to specify vector_direction_of_unit_cut?
vector_direction_of_unit_cut = [];
if nargin>=4
    temp = varargin{3};
    if ~isempty(temp)
        vector_direction_of_unit_cut = temp;
    end
end

% Does user want to specify flag_vertexIsNonConvex?
flag_vertexIsNonConvex = [];
if nargin>=5
    temp = varargin{4};
    if ~isempty(temp)
        flag_vertexIsNonConvex = temp;
    end
end

% Does user want to specify flag_plotEdgeGhostlines?
flag_plotEdgeGhostlines = [];
if nargin>=6
    temp = varargin{5};
    if ~isempty(temp)
        flag_plotEdgeGhostlines = temp;
    end
end

% Does user want to specify flag_plotVertexProjectionGhostlines?
flag_plotVertexProjectionGhostlines = [];
if nargin>=7
    temp = varargin{6};
    if ~isempty(temp)
        flag_plotVertexProjectionGhostlines = temp;
    end
end


% Fill in defaults
plot_formatting.vertices_plot.style = '.-';
plot_formatting.vertices_plot.LineWidth = 2;
plot_formatting.vertices_plot.MarkerSize = 20;
plot_formatting.vertices_plot.Color = [0 0 1];
plot_formatting.vertices_plot.vertexLabelsColor = [0 1 0];
plot_formatting.vertices_plot.edgeLabelsColor = [0 0 1];

plot_formatting.edgeGhostLines_plot.style = '-';
plot_formatting.edgeGhostLines_plot.LineWidth = 0.5;
plot_formatting.edgeGhostLines_plot.MarkerSize = 0.1;
plot_formatting.edgeGhostLines_plot.Color = 0.7*[1 1 1];

plot_formatting.vertexProjectionGhostLines_plot.style = '-';
plot_formatting.vertexProjectionGhostLines_plot.LineWidth = 0.5;
plot_formatting.vertexProjectionGhostLines_plot.MarkerSize = 0.1;
plot_formatting.vertexProjectionGhostLines_plot.Color = 0.7*[0 1 0];

plot_formatting.unitNormalVectors_plot.style = 'r';
plot_formatting.unitNormalVectors_plot.LineWidth = 0.5;
plot_formatting.unitNormalVectors_plot.MarkerSize = 0.1;
plot_formatting.unitNormalVectors_plot.Color = [1 0 0];

plot_formatting.unitVertexProjectionVectors_plot.style = 'g';
plot_formatting.unitVertexProjectionVectors_plot.LineWidth = 3;
plot_formatting.unitVertexProjectionVectors_plot.MarkerSize = 0.1;
plot_formatting.unitVertexProjectionVectors_plot.Color = [0 1 0];

plot_formatting.vectorDirectionOfUnitCut_plot.style = '-';
plot_formatting.vectorDirectionOfUnitCut_plot.LineWidth = 2;
plot_formatting.vectorDirectionOfUnitCut_plot.MarkerSize = 0.1;
plot_formatting.vectorDirectionOfUnitCut_plot.Color = [0 0.5 0];

plot_formatting.vertexIsNonConvex_plot.style = '.';
plot_formatting.vertexIsNonConvex_plot.LineWidth = 2;
plot_formatting.vertexIsNonConvex_plot.MarkerSize = 20;
plot_formatting.vertexIsNonConvex_plot.Color = [1 0 0];

% Does user want to specify plot_formatting?
if nargin>=8
    temp = varargin{7};
    if ~isempty(temp)
        plot_formatting = fcn_INTERNAL_copyStructIntoStruct(plot_formatting,temp);
    end
end

% Does user want to show the plots?
flag_do_plot = 1; % Default is ALWAYS plotting
if  9 == nargin 
    temp = varargin{end}; % Last argument is always figure number
    if ~isempty(temp) % Make sure the user is not giving empty input
        fig_num = temp;
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

% Preliminary calculations go here, before plotting
NumUniqueVerticies = length(vertices(:,1));

% Find size of vertex domain
max_XY = max(vertices);
min_XY = min(vertices);
sizePlot = max(max_XY) - min(min_XY);
nudge = sizePlot*0.006;

% Find the modpoints for each vertex
midpoints = (vertices(2:end,:)+vertices(1:end-1,:))/2;
midpoints = [midpoints; midpoints(1,:)]; % Repeat first row, to last, to match how point is similarly repeated

[INTERNAL_unit_normal_vectors, INTERNAL_unit_vertex_projection_vectors]  = ...
    fcn_VSkel_polytopeFindUnitDirectionVectors(vertices,-1);

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
    h_fig = figure(fig_num);
    flag_rescale_axis = 0;
    if isempty(get(h_fig,'Children'))
        flag_rescale_axis = 1;
    end    

    % grid on
    % grid minor
    hold on
    axis equal   

    % Find size of vertex domain
    if flag_rescale_axis
        axis_range_x = max_XY(1,1)-min_XY(1,1);
        axis_range_y = max_XY(1,2)-min_XY(1,2);
        percent_larger = 0.3;
        axis([min_XY(1,1)-percent_larger*axis_range_x, max_XY(1,1)+percent_larger*axis_range_x,  min_XY(1,2)-percent_larger*axis_range_y, max_XY(1,2)+percent_larger*axis_range_y]);    
    end
    goodAxis = axis;
    

    % Plot the polytope in dots connected by lines
    plot(vertices(:,1),vertices(:,2),plot_formatting.vertices_plot.style,'Linewidth',plot_formatting.vertices_plot.LineWidth, 'MarkerSize',plot_formatting.vertices_plot.MarkerSize, 'Color',plot_formatting.vertices_plot.Color);

    % Label the vertices with their numbers
    for ith_vertex = 1:length(vertices(:,1))-1
        text(vertices(ith_vertex,1)+nudge,vertices(ith_vertex,2),...
            sprintf('%.0d',ith_vertex),'Color',plot_formatting.vertices_plot.vertexLabelsColor);
    end

    % Label the edges with their numbers
    for ith_edge = 1:length(vertices(:,1))-1
        text(midpoints(ith_edge,1)+nudge, midpoints(ith_edge,2),...
            sprintf('%.0d',ith_edge),'Color',plot_formatting.vertices_plot.edgeLabelsColor);
    end

    % Plot the edge "ghostlines"
    if ~isempty(flag_plotEdgeGhostlines) && (1==flag_plotEdgeGhostlines)
        unit_tangent_vectors = INTERNAL_unit_normal_vectors*[0 1; -1 0];
        for ith_vertex = 1:length(vertices(:,1))-1
            ghostEnds = [...
                vertices(ith_vertex,:)+2*sizePlot*unit_tangent_vectors(ith_vertex,:);
                vertices(ith_vertex,:)-2*sizePlot*unit_tangent_vectors(ith_vertex,:);
                ];
            plot(ghostEnds(:,1),ghostEnds(:,2),plot_formatting.edgeGhostLines_plot.style,'Linewidth',plot_formatting.edgeGhostLines_plot.LineWidth, 'MarkerSize',plot_formatting.edgeGhostLines_plot.MarkerSize, 'Color',plot_formatting.edgeGhostLines_plot.Color);

        end
    end

    % Plot the vertex "ghostlines"
    if ~isempty(flag_plotVertexProjectionGhostlines)  && (1==flag_plotVertexProjectionGhostlines)
        for ith_vertex = 1:length(vertices(:,1))-1
            ghostEnds = [...
                vertices(ith_vertex,:)+0*sizePlot*INTERNAL_unit_vertex_projection_vectors(ith_vertex,:);
                vertices(ith_vertex,:)+2*sizePlot*INTERNAL_unit_vertex_projection_vectors(ith_vertex,:);
                ];
            plot(ghostEnds(:,1),ghostEnds(:,2),plot_formatting.vertexProjectionGhostLines_plot.style,'Linewidth',plot_formatting.vertexProjectionGhostLines_plot.LineWidth, 'MarkerSize',plot_formatting.vertexProjectionGhostLines_plot.MarkerSize, 'Color',plot_formatting.vertexProjectionGhostLines_plot.Color);

        end
    end

    % Draw the unit vectors
    if ~isempty(unit_normal_vectors)
        quiver(midpoints(1:end-1,1),midpoints(1:end-1,2),unit_normal_vectors(1:end-1,1),unit_normal_vectors(1:end-1,2),0, ...
            plot_formatting.unitNormalVectors_plot.style,'Linewidth',plot_formatting.unitNormalVectors_plot.LineWidth, 'MarkerSize',plot_formatting.unitNormalVectors_plot.MarkerSize, 'Color',plot_formatting.unitNormalVectors_plot.Color);
    end

    % Draw the vertex_projection_vectors 
    if ~isempty(unit_vertex_projection_vectors)
        quiver(vertices(1:end-1,1),vertices(1:end-1,2), unit_vertex_projection_vectors(1:end-1,1),unit_vertex_projection_vectors(1:end-1,2),0, ...
            plot_formatting.unitVertexProjectionVectors_plot.style,'Linewidth',plot_formatting.unitVertexProjectionVectors_plot.LineWidth, 'MarkerSize',plot_formatting.unitVertexProjectionVectors_plot.MarkerSize, 'Color',plot_formatting.unitVertexProjectionVectors_plot.Color);
    end

    % Draw the vector_direction_of_unit_cut
    if ~isempty(vector_direction_of_unit_cut)
        quiver(vertices(1:end-1,1),vertices(1:end-1,2), vector_direction_of_unit_cut(1:end-1,1),vector_direction_of_unit_cut(1:end-1,2),0,...
            plot_formatting.vectorDirectionOfUnitCut_plot.style,'Linewidth',plot_formatting.vectorDirectionOfUnitCut_plot.LineWidth, 'MarkerSize',plot_formatting.vectorDirectionOfUnitCut_plot.MarkerSize, 'Color',plot_formatting.vectorDirectionOfUnitCut_plot.Color);
    end

    % Label any non-convex verticies
    if ~isempty(flag_vertexIsNonConvex)
        bad_verticies = find(flag_vertexIsNonConvex);
        plot(vertices(bad_verticies,1),vertices(bad_verticies,2),...
            plot_formatting.vertexIsNonConvex_plot.style,'Linewidth',plot_formatting.vertexIsNonConvex_plot.LineWidth, 'MarkerSize',plot_formatting.vertexIsNonConvex_plot.MarkerSize, 'Color',plot_formatting.vertexIsNonConvex_plot.Color);
    end

    axis(goodAxis);



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
%% fcn_INTERNAL_copyStructIntoStruct
function output_struct = fcn_INTERNAL_copyStructIntoStruct(input_struct,temp)
output_struct = input_struct;
if isstruct(temp)
    names_to_copy = fieldnames(temp);
    for ith_field = 1:length(names_to_copy)
        fieldname = names_to_copy{ith_field};
        output_struct.(fieldname) = fcn_INTERNAL_copyStructIntoStruct(input_struct.(fieldname),temp.(fieldname));
    end
else
    output_struct = temp;
end

end % Ends fcn_INTERNAL_copyStructIntoStruct