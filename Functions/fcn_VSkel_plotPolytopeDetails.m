function h_fig = fcn_VSkel_plotPolytopeDetails(vertices, varargin)

%% fcn_VSkel_plotPolytopeDetails
% plots a single poloytpe allowing added user-defined details
%
% FORMAT:
%
% h_fig =  fcn_VSkel_plotPolytopeDetails(...
%        polytopeStructure,...
%        (faceVectorsToPlot), ...  % faceVectorsToPlot
%        (vertexVectorsToPlot), ...  % vertexVectorsToPlot
%        (plot_formatting),... % plot_formatting
%        (fig_num));  % fig_num
%
% INPUTS:
%
%     vertices: a (M+1)-by-2 matrix of xy points with row1 = rowm+1, where
%     M is the number of the individual polytope vertices. If a
%     polytopeStructure type is given for vertices, it is plotted using
%     polytopeStructure formatting (see polytopeFillStructureFromVertices
%     for details
%
%     (OPTIONAL INPUTS)
%
%     faceVectorsToPlot: an F-by-2 or F-by-3 matrix of user-defined vectors
%     that are anchored to the midpoint of each face, where F is the number
%     of faces and 2 or 3 refer to either 2D or 3D points.
%
%     vertexVectorsToPlot: a V-by-2 or V-by-3 matrix of user-defined
%     vectors that are anchored at each vertex, where V is the number of
%     verticies and 2 or 3 refer to either 2D or 3D points.
%
%     plot_formatting: a structure specifying the plot style. For any
%     fields left empty, defaults are used. The defaults are
%
%         plot_formatting.vertices_plot.style = '.-';
%         plot_formatting.vertices_plot.LineWidth = 2;
%         plot_formatting.vertices_plot.MarkerSize = 20;
%         plot_formatting.vertices_plot.Color = [0 0 1];
%         plot_formatting.vertices_plot.vertexLabels_Color = 0*[1 1 1];
%         plot_formatting.vertices_plot.faceLabels_Color = [0 0 1];
%         
%         plot_formatting.faceVectorsToPlot_plot.style = '-';
%         plot_formatting.faceVectorsToPlot_plot.LineWidth = 0.5;
%         plot_formatting.faceVectorsToPlot_plot.MarkerSize = 0.1;
%         plot_formatting.faceVectorsToPlot_plot.Color = [1 0 0];
%         
%         plot_formatting.vertexVectorsToPlot_plot.style = '-';
%         plot_formatting.vertexVectorsToPlot_plot.LineWidth = 3;
%         plot_formatting.vertexVectorsToPlot_plot.MarkerSize = 0.1;
%         plot_formatting.vertexVectorsToPlot_plot.Color = [0 1 0];
%         
%         plot_formatting.edgeGhostLines_flagOn = 0; % Flag to indicate that edgeGhostLines should be plotted
%         plot_formatting.edgeGhostLines_plot.style = '-';
%         plot_formatting.edgeGhostLines_plot.LineWidth = 0.5;
%         plot_formatting.edgeGhostLines_plot.MarkerSize = 0.1;
%         plot_formatting.edgeGhostLines_plot.Color = 0.7*[1 1 1];
%         
%         plot_formatting.vertexProjectionGhostLines_flagOn = 0; % Flag to indicate that vertexProjectionGhostLines should be plotted
%         plot_formatting.vertexProjectionGhostLines_plot.style = '-';
%         plot_formatting.vertexProjectionGhostLines_plot.LineWidth = 0.5;
%         plot_formatting.vertexProjectionGhostLines_plot.MarkerSize = 0.1;
%         plot_formatting.vertexProjectionGhostLines_plot.Color = 0.7*[0 1 0];
%         
%         plot_formatting.vertexIsNonConvex_flagOn = 0; % Flag to indicate that non-convex points should be plotted
%         plot_formatting.vertexIsNonConvex_plot.Marker  = '.';
%         plot_formatting.vertexIsNonConvex_plot.LineWidth = 2;
%         plot_formatting.vertexIsNonConvex_plot.MarkerSize = 20;
%         plot_formatting.vertexIsNonConvex_plot.Color = [1 0 0];
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
% 2025_05_15 by Sean Brennan
% -- added case where vertices can be only one point
% 2025_05_27 by Sean Brennan
% -- merged flags and plot details for simplicity


% TO DO
% -- none

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
flag_max_speed = 0;
if (nargin==5 && isequal(varargin{end},-1))
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
        narginchk(1,5);

        % Check the vertices input
        if ~isstruct(vertices)
            fcn_DebugTools_checkInputsToFunctions(...
                vertices, '2or3column_of_numbers');
        end
    end
end

% Is the user entering a polytopeStructure?
if isstruct(vertices)
    polytopeStructure = vertices;
else
    polytopeStructure = [];
end

% Does user want to specify faceVectorsToPlot?
faceVectorsToPlot = [];
if nargin>=2
    temp = varargin{1};
    if ~isempty(temp)
        faceVectorsToPlot = temp;
    end
end


% Does user want to specify vertexVectorsToPlot?
vertexVectorsToPlot = [];
if nargin>=3
    temp = varargin{2};
    if ~isempty(temp)
        vertexVectorsToPlot = temp;
    end
end



% Fill in defaults
plot_formatting.vertices_plot.style = '.-';
plot_formatting.vertices_plot.LineWidth = 2;
plot_formatting.vertices_plot.MarkerSize = 20;
plot_formatting.vertices_plot.Color = [0 0 1];
plot_formatting.vertices_plot.vertexLabels_flagOn = 0; % Flag to indicate that vertexLabels should be shown
plot_formatting.vertices_plot.edgeLabels_flagOn = 0; % Flag to indicate that edgeLabels should be shown
plot_formatting.vertices_plot.vertexLabels_Color = 0*[1 1 1];
plot_formatting.vertices_plot.edgeLabels_Color = 0*[1 0 0];
plot_formatting.vertices_plot.faceLabels_Color = [0 0 1];
plot_formatting.vertices_plot.faceLabels_flagOn = 0; % Flag to indicate that faceLabels should be shown

plot_formatting.faceVectorsToPlot_plot.style = '-';
plot_formatting.faceVectorsToPlot_plot.LineWidth = 0.5;
plot_formatting.faceVectorsToPlot_plot.MarkerSize = 0.1;
plot_formatting.faceVectorsToPlot_plot.Color = [1 0 0];

plot_formatting.vertexVectorsToPlot_plot.style = '-';
plot_formatting.vertexVectorsToPlot_plot.LineWidth = 3;
plot_formatting.vertexVectorsToPlot_plot.MarkerSize = 0.1;
plot_formatting.vertexVectorsToPlot_plot.Color = [0 1 0];

plot_formatting.edgeGhostLines_flagOn = 0; % Flag to indicate that edgeGhostLines should be plotted
plot_formatting.edgeGhostLines_plot.style = '-';
plot_formatting.edgeGhostLines_plot.LineWidth = 0.5;
plot_formatting.edgeGhostLines_plot.MarkerSize = 0.1;
plot_formatting.edgeGhostLines_plot.Color = 0.7*[1 1 1];

plot_formatting.vertexProjectionGhostLines_flagOn = 0; % Flag to indicate that vertexProjectionGhostLines should be plotted
plot_formatting.vertexProjectionGhostLines_plot.style = '-';
plot_formatting.vertexProjectionGhostLines_plot.LineWidth = 0.5;
plot_formatting.vertexProjectionGhostLines_plot.MarkerSize = 0.1;
plot_formatting.vertexProjectionGhostLines_plot.Color = 0.7*[0 1 0];

plot_formatting.vertexIsNonConvex_flagOn = 0; % Flag to indicate that non-convex points should be plotted
plot_formatting.vertexIsNonConvex_plot.Marker  = '.';
plot_formatting.vertexIsNonConvex_plot.LineWidth = 2;
plot_formatting.vertexIsNonConvex_plot.MarkerSize = 20;
plot_formatting.vertexIsNonConvex_plot.Color = [1 0 0];

% Does user want to specify plot_formatting?
if nargin>=4
    temp = varargin{3};
    if ~isempty(temp)
        plot_formatting = fcn_INTERNAL_copyStructIntoStruct(plot_formatting,temp);
    end
end

% Does user want to show the plots?
flag_do_plot = 1; % Default is ALWAYS plotting
if  5 == nargin
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

if ~isempty(polytopeStructure)
    vertices = polytopeStructure.polyPatch.Vertices;
end

% Is this 2D or 3D?
dimension_of_points = length(vertices(1,:));

% How many things are there?
Npolytopes = length(polytopeStructure.subPolyPatch); %#ok<NASGU>
Nvertices  = length(vertices(:,1));
Nfaces = length(polytopeStructure.polyPatch.Faces(:,1));

% Find size of vertex domain
max_vertexValues = max(vertices);
min_vertexValues = min(vertices);
sizePlot = max(max_vertexValues) - min(min_vertexValues);
nudge = sizePlot*0.006;

% % Find the modpoints for each Edge
% if 3==dimension_of_points
%     midpointsEdges = zeros(Nfaces,dimension_of_points);
%     Nmidpoints = Nfaces;
%     facesOrEdges = polytopeStructure.polyPatch.Faces;
% else
%     midpointsFaces = zeros(Nvertices,dimension_of_points);
%     Nmidpoints = Nvertices;
%     facesOrEdges = [];
%     for ith_face = 1:length(polytopeStructure.subPolyPatch)
%         facesOrEdges = [facesOrEdges; polytopeStructure.subPolyPatch(ith_face).Faces]; %#ok<AGROW>
%     end
% end


% Find the modpoints for each Face
if 3==dimension_of_points
    midpointsFaces = zeros(Nfaces,dimension_of_points);
    Nmidpoints = Nfaces;
    facesOrEdges = polytopeStructure.polyPatch.Faces;
else
    midpointsFaces = zeros(Nvertices,dimension_of_points);
    Nmidpoints = Nvertices;
    facesOrEdges = [];
    for ith_face = 1:length(polytopeStructure.subPolyPatch)
        facesOrEdges = [facesOrEdges; polytopeStructure.subPolyPatch(ith_face).Faces]; %#ok<AGROW>
    end
end

for ith_face = 1:Nmidpoints
    indicesPointsInFace   = (facesOrEdges(ith_face,:))';
    realIndices           = indicesPointsInFace(~isnan(indicesPointsInFace));
    pointsInFace          = vertices(realIndices,:);
    midpointsFaces(ith_face,:) = mean(pointsInFace,1,'omitmissing');
end

% Find the projection vectors?
if plot_formatting.edgeGhostLines_flagOn == 1 || plot_formatting.vertexProjectionGhostLines_flagOn == 1  ||  plot_formatting.vertexIsNonConvex_flagOn == 1 
    [INTERNAL_unit_normal_vectors, INTERNAL_vector_direction_of_unit_cut, INTERNAL_flag_vertexIsNonConvex] = ...
        fcn_VSkel_polytopeFindUnitDirectionVectors(polytopeStructure,-1);

    if ~iscell(INTERNAL_vector_direction_of_unit_cut)
        lengths = sum(INTERNAL_vector_direction_of_unit_cut.^2,2).^0.5;
        INTERNAL_unit_vertex_vectors = INTERNAL_vector_direction_of_unit_cut./lengths;
    end

end

%% Plot results?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____       _
%  |  __ \     | |
%  | |  | | ___| |__  _   _  _
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
    if isempty(get(h_fig,'Children')) || strcmp(get(get(h_fig,'Children'),'TileArrangement'),'flow')
        flag_rescale_axis = 1; 
    end

    % grid on
    % grid minor
    hold on


    % Find size of vertex domain
    if flag_rescale_axis
        percent_larger = 0.3;
        axis_range = max_vertexValues - min_vertexValues;
        if (0==axis_range(1,1))
            axis_range(1,1) = 2/percent_larger;
        end
        if (0==axis_range(1,2))
            axis_range(1,2) = 2/percent_larger;
        end
        if dimension_of_points==3 && (0==axis_range(1,3))
            axis_range(1,3) = 2/percent_larger;
        end

        
        % Force the axis to be equal
        min_vertexValuesInPlot = min(min_vertexValues);
        max_vertexValuesInPlot = max(max_vertexValues);

        % Stretch the axes
        stretched_min_vertexValues = min_vertexValuesInPlot - percent_larger.*axis_range;
        stretched_max_vertexValues = max_vertexValuesInPlot + percent_larger.*axis_range;
        axesTogether = [stretched_min_vertexValues; stretched_max_vertexValues];
        newAxis = reshape(axesTogether, 1, []);
        axis(newAxis);

    end
    goodAxis = axis;


    if ~isempty(polytopeStructure)
        patch(polytopeStructure.polyPatch);

        % % Plot all the polytopes with shading
        % for ith_polytope = 1:Npolytopes
        %     patch(polytopeStructure.subPolyPatch(ith_polytope));
        % end

        xlabel('X');
        ylabel('Y');
        zlabel('Z');

    else

        if dimension_of_points==3 
            error('Need to fix this.');
        end
        % Plot the polytope in dots connected by lines
        fcn_INTERNAL_plotND(vertices,...
            'LineStyle', plot_formatting.vertices_plot.style,'Linewidth',plot_formatting.vertices_plot.LineWidth, 'MarkerSize',plot_formatting.vertices_plot.MarkerSize, 'Color',plot_formatting.vertices_plot.Color);
    end

    % Label the vertices with their numbers?
    if plot_formatting.vertices_plot.vertexLabels_flagOn == 1
        nudgedVertices = vertices;
        nudgedVertices(:,1) = nudgedVertices(:,1)+nudge;
        for ith_vertex = 1:Nvertices
            fcn_INTERNAL_textND(nudgedVertices(ith_vertex,:),...
                sprintf('%.0d',ith_vertex),'Color',plot_formatting.vertices_plot.vertexLabels_Color);

        end
    end

    % % Label the edges with their numbers?
    % if plot_formatting.vertices_plot.edgeLabels_flagOn == 1
    %     nudgedEdges = vertices;
    %     nudgedEdges(:,1) = nudgedEdges(:,1)+nudge;
    %     for ith_vertex = 1:Nvertices
    %         fcn_INTERNAL_textND(nudgedVertices(ith_vertex,:),...
    %             sprintf('%.0d',ith_vertex),'Color',plot_formatting.vertices_plot.vertexLabels_Color);
    % 
    %     end
    % end


    % Label the faces with their numbers?    
    if plot_formatting.vertices_plot.faceLabels_flagOn == 1
        nudgedMidpoints = midpointsFaces;
        nudgedMidpoints(:,1) = nudgedMidpoints(:,1) + nudge;
        for ith_edge = 1:Nmidpoints
            fcn_INTERNAL_textND(nudgedMidpoints(ith_edge,:),...
                sprintf('%.0d',ith_edge),'Color',plot_formatting.vertices_plot.faceLabels_Color);
        end
    end

    % Plot the edge "ghostlines"?
    if plot_formatting.edgeGhostLines_flagOn == 1
        unit_tangent_vectors = INTERNAL_unit_normal_vectors*[0 1; -1 0];
        for ith_vertex = 1:length(vertices(:,1))-1
            ghostEnds = [...
                vertices(ith_vertex,:)+2*sizePlot*unit_tangent_vectors(ith_vertex,:);
                vertices(ith_vertex,:)-2*sizePlot*unit_tangent_vectors(ith_vertex,:);
                ];
            fcn_INTERNAL_plotND(ghostEnds, ...
                'LineStyle',plot_formatting.edgeGhostLines_plot.style,'Linewidth',plot_formatting.edgeGhostLines_plot.LineWidth, 'MarkerSize',plot_formatting.edgeGhostLines_plot.MarkerSize, 'Color',plot_formatting.edgeGhostLines_plot.Color);

        end
    end

    % Plot the vertex "ghostlines"?
    if  plot_formatting.vertexProjectionGhostLines_flagOn == 1
        for ith_vertex = 1:length(vertices(:,1))-1
            ghostEnds = [...
                vertices(ith_vertex,:)+0*sizePlot*INTERNAL_unit_vertex_vectors(ith_vertex,:);
                vertices(ith_vertex,:)+2*sizePlot*INTERNAL_unit_vertex_vectors(ith_vertex,:);
                ];
            fcn_INTERNAL_plotND(ghostEnds, ...
                'LineStyle', plot_formatting.vertexProjectionGhostLines_plot.style,'Linewidth',plot_formatting.vertexProjectionGhostLines_plot.LineWidth, 'MarkerSize',plot_formatting.vertexProjectionGhostLines_plot.MarkerSize, 'Color',plot_formatting.vertexProjectionGhostLines_plot.Color);

        end
    end

    % Draw the face vectors?
    if ~isempty(faceVectorsToPlot)
        fcn_INTERNAL_quiverND(midpointsFaces,faceVectorsToPlot, ...
            'LineStyle', plot_formatting.faceVectorsToPlot_plot.style,'Linewidth',plot_formatting.faceVectorsToPlot_plot.LineWidth, 'MarkerSize',plot_formatting.faceVectorsToPlot_plot.MarkerSize, 'Color',plot_formatting.faceVectorsToPlot_plot.Color);
    end

    % Draw the vertex vectors?
    if ~isempty(vertexVectorsToPlot)
        fcn_INTERNAL_quiverND(vertices, vertexVectorsToPlot, ...
            'LineStyle', plot_formatting.vertexVectorsToPlot_plot.style,'Linewidth',plot_formatting.vertexVectorsToPlot_plot.LineWidth, 'MarkerSize',plot_formatting.vertexVectorsToPlot_plot.MarkerSize, 'Color',plot_formatting.vertexVectorsToPlot_plot.Color);
    end

    % Label any non-convex verticies?
    if plot_formatting.vertexIsNonConvex_flagOn == 1 
        bad_verticies = INTERNAL_flag_vertexIsNonConvex;
        if any(bad_verticies~=0)
            fcn_INTERNAL_plotND(vertices(bad_verticies,:),...
                'Marker', plot_formatting.vertexIsNonConvex_plot.Marker ,'Linewidth',plot_formatting.vertexIsNonConvex_plot.LineWidth, 'MarkerSize',plot_formatting.vertexIsNonConvex_plot.MarkerSize, 'Color',plot_formatting.vertexIsNonConvex_plot.Color);
        end
    end

    axis(goodAxis);
    axis square;

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


%% fcn_INTERNAL_plotND
function fcn_INTERNAL_plotND(vertices, varargin)

% Is this 2D or 3D?
dimension_of_points = length(vertices(1,:));
if 2==dimension_of_points
    h_plot = plot(vertices(:,1), vertices(:,2));
else
    h_plot = plot3(vertices(:,1), vertices(:,2));
end
assert(mod(length(varargin),2)==0);
Narguments = length(varargin)/2;
for ith_argument = 1:Narguments
    first_argument_index = (ith_argument-1)*2 + 1;
    second_argument_index = (ith_argument-1)*2 + 2;
    field_to_set = varargin{first_argument_index};
    setting = varargin{second_argument_index};
    set(h_plot,field_to_set,setting);
end


end % End fcn_INTERNAL_plotNDS

%% fcn_INTERNAL_quiverND
function fcn_INTERNAL_quiverND(vertices, vectors, varargin)

% Is this 2D or 3D?
dimension_of_points = length(vertices(1,:));
if 2==dimension_of_points
    h_quiver = quiver(vertices(:,1),vertices(:,2), vectors(:,1), vectors(:,2), 0);
else
    if ~iscell(vectors)
        h_quiver = quiver3(vertices(:,1),vertices(:,2), vertices(:,3), vectors(:,1), vectors(:,2), vectors(:,3), 0);
        assert(mod(length(varargin),2)==0);
        Narguments = length(varargin)/2;
        for ith_argument = 1:Narguments
            first_argument_index = (ith_argument-1)*2 + 1;
            second_argument_index = (ith_argument-1)*2 + 2;
            field_to_set = varargin{first_argument_index};
            setting = varargin{second_argument_index};
            set(h_quiver,field_to_set,setting);
        end
    else
        for ith_vertex = 1:length(vertices(:,1))
            theseVectors = vectors{ith_vertex};
            for ith_vector = 1:length(theseVectors(:,1))
                h_quiver = quiver3(vertices(ith_vertex,1),vertices(ith_vertex,2), vertices(ith_vertex,3), theseVectors(ith_vector,1), theseVectors(ith_vector,2), theseVectors(ith_vector,3), 0);
                assert(mod(length(varargin),2)==0);
                Narguments = length(varargin)/2;
                for ith_argument = 1:Narguments
                    first_argument_index = (ith_argument-1)*2 + 1;
                    second_argument_index = (ith_argument-1)*2 + 2;
                    field_to_set = varargin{first_argument_index};
                    setting = varargin{second_argument_index};
                    set(h_quiver,field_to_set,setting);
                end
            end
        end

    end        
end

end % End fcn_INTERNAL_quiverND


%% fcn_INTERNAL_textND
function fcn_INTERNAL_textND(vertices, string, varargin)

% Is this 2D or 3D?
dimension_of_points = length(vertices(1,:));
if 2==dimension_of_points
    h_text = text(vertices(:,1),vertices(:,2), string);
else
    h_text = text(vertices(:,1),vertices(:,2), vertices(:,3), string);    
end

assert(mod(length(varargin),2)==0);
Narguments = length(varargin)/2;
for ith_argument = 1:Narguments
    first_argument_index = (ith_argument-1)*2 + 1;
    second_argument_index = (ith_argument-1)*2 + 2;
    field_to_set = varargin{first_argument_index};
    setting = varargin{second_argument_index};
    set(h_text,field_to_set,setting);
end

end % End fcn_INTERNAL_quiverND

