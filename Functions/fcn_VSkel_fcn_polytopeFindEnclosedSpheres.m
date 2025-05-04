function [sphereRadii, ...
    definingEdges] = ...
    fcn_VSkel_fcn_polytopeFindEnclosedSpheres(vertices, unit_normal_vectors, vertex_projection_vectors, varargin)

%% fcn_VSkel_fcn_polytopeFindEnclosedSpheres
% finds the enclosed sphere radii for each vertex point, to each other edge
%
% FORMAT:
%
% [sphereRadii, definingEdges] = ...
%     fcn_VSkel_fcn_polytopeFindEnclosedSpheres(vertices, unit_normal_vectors, vertex_projection_vectors, (fig_num))
%
% INPUTS:
%
%     vertices: a (M+1)-by-2 matrix of xy points with row1 = rowm+1, where
%         M is the number of the individual polytope vertices
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
%     sphereRadii: a cell array of dimension N containing, in each cell, an
%     array of radii. In each cell array, the are Mx1 radii vectors, where
%     M is = N-2 and N is the number of vertices. The radii are ordered so
%     that the first radii cell array corresponds to the first vertex, etc.
%
%     definingEdges: a cell array of dimension N containing, in
%     each cell, an array of which edges define each vertex. In each cell
%     array, there are Mx1 defining edges, where M is = N-2 and N is the
%     number of vertices. The defining edges match the radii ordering, e.g.
%     vertex 2's 3rd sphereRadii edge interaction ID will be in cell array
%     2, in the 3rd row.
%
% DEPENDENCIES:
%
%     fcn_DebugTools_checkInputsToFunctions
%
% % EXAMPLES:
% For additional examples, see: script_test_fcn_VSkel_fcn_polytopeFindEnclosedSpheres
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
% -- speed up by checking off which edges were checked previously, and not
% repeating the calculations again for these. For example, in a 4-sided 2D
% polytope, vertex 1 contains edges 41. If this hits edge 3, the result
% will be the same as vertex 4 which is edges 34 hitting edge 1. They all
% contain edges 134.

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

if flag_check_inputs
    % Are there the right number of inputs?
    if nargin < 3 || nargin > 4
        error('Incorrect number of input arguments')
    end

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

    % Check the vertex_projection_vectors input
    fcn_DebugTools_checkInputsToFunctions(...
        vertex_projection_vectors, '2or3column_of_numbers',NumUniqueVerticies);

end


% Does user want to show the plots?
flag_do_plot = 0; % Default is no plotting
if  4 == nargin && (0==flag_max_speed) % Only create a figure if NOT maximizing speed
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
    sphereRadii = cell(NumUniqueVerticies,1);
    sphereCenters = cell(NumUniqueVerticies,1);
    definingEdges = cell(NumUniqueVerticies,1);

    % For each vertex, solve for the radii to all the non-participating
    % edges. Non-participating edges are those that are not to either side
    % of this vertex.
    for ith_vertex = 1:NumUniqueVerticies
        % Find needed values for current vertex
        vertexProjection = vertex_projection_vectors(ith_vertex,:);
        vertexNormal = unit_normal_vectors(ith_vertex,:);
        vertexPoint  = vertices(ith_vertex,:);

        % Find need values for previous vertex
        if ith_vertex == 1
            previousVertexIndex = NumUniqueVerticies;
        else
            previousVertexIndex = ith_edge-1;
        end
        previousVertexNormal = unit_normal_vectors(previousVertexIndex,:);
        

        % Initialize arrays
        radiiFromVertexToEdge  = zeros(NumUniqueVerticies-2,1);
        edgesConstrainingRadii = zeros(NumUniqueVerticies-2,1);
        sphereCenterArray      = zeros(NumUniqueVerticies-2,2);

        % The edges tested start with the edge ahead of this vertex's edge (which is the one
        % to the right of the vertex), and move up until the one before
        % this vertex, to the left. This will omit 2 verticies. Including
        % the 2 omitted, and the one we are on, we go up by
        % NumUniqueVerticies - 3. We convert number back into real edge
        % numbering. For example, if there are 4 edges, edge 5 is really
        % edge 1, edge 6 is really edge 2, etc. The mod function winds the
        % for-loop numbering back into real edge numbering.
        edgesToTest = (ith_vertex:(ith_vertex+NumUniqueVerticies-3))';
        edgesToTest = mod(edgesToTest,NumUniqueVerticies)+1;

        % Loop through non-participating edges.
        for ith_edge = 1:length(edgesToTest)
            current_edge = edgesToTest(ith_edge,1);
            edgePointStart = vertices(current_edge,:);
            edgeNormal = unit_normal_vectors(current_edge,:);



            % See derivation in PPT documentation for the equations below.
            % Requires diagrams to explain where they come from
            d = sum((vertexPoint - edgePointStart).*(edgeNormal),2) / sum(((vertexNormal - edgeNormal).*vertexProjection),2);
            r = d*sum(vertexProjection.*vertexNormal,2);

            % The sphere's center will be the projection from the vertex
            % point by distance d
            sphereCenter = vertexPoint + vertexProjection*d;

            if 1==1
                % Check to see if point is on edge by projecting it back toward
                % the edge by the radius times the edgeNormal
                testPoint = sphereCenter - edgeNormal*r;
                edgePointEnd = vertices(current_edge+1,:);

                % For debugging
                if 1==flag_do_debug
                    figure(6666);
                    clf;
                    % Plot the polytope in black dots connected by lines
                    plot(vertices(:,1),vertices(:,2),'b.-','Linewidth',2, 'MarkerSize',10);
                    hold on;


                    % Plot the circle center as a large dot, and store the color
                    h_fig = plot(sphereCenter(1,1),sphereCenter(1,2),'.','MarkerSize',20);
                    colorUsed = get(h_fig,'Color');

                    % Plot the circle boundary in same color
                    fcn_geometry_plotCircle(sphereCenter,r,colorUsed,6666);

                    % Plot the edge being tested in same color
                    edgeData = [edgePointStart; edgePointEnd];
                    plot(edgeData(:,1),edgeData(:,2),'-','Linewidth',2, 'Color', colorUsed);

                    % Plot the current vertex in same color
                    plot(vertices(ith_vertex,1),vertices(ith_vertex,2),'.','MarkerSize',50, 'Color', colorUsed);


                    plot(testPoint(1,1),testPoint(1,2),'.','MarkerSize',50)
                end


                % Check if test point is on edge
                isOnEdge = fcn_INTERNAL_isPointOnEdge(edgePointStart,edgePointEnd,testPoint);

                % Check if test point is in positive direction for both
                % base orthogonal projections
                isWithin = fcn_INTERNAL_isPointWithinVertexProjections(...
                    testPoint, vertexPoint,  vertexProjection,...
                    vertexNormal, previousVertexNormal);

                if 1==flag_do_debug
                    figure(6666);
                    title(sprintf('Vertex: %.0d, isOnEdge: %.0f, isWithin: %.0f',ith_vertex, isOnEdge, isWithin));
                end

                if ~isOnEdge && ~isWithin
                    r = nan;
                    sphereCenter = [nan nan];
                end
            end
            % Save results to matricies
            radiiFromVertexToEdge(ith_edge,1)  = r;
            sphereCenterArray(ith_edge,:) = sphereCenter;
            edgesConstrainingRadii(ith_edge,1) = current_edge;
        end

        % Save matricies to cell arrays for this vertex
        sphereRadii{ith_vertex}   = radiiFromVertexToEdge;
        sphereCenters{ith_vertex} = sphereCenterArray;
        definingEdges{ith_vertex} = edgesConstrainingRadii;
    end

    % Repeat last value as first, since this is the repeated 1st vertex
    sphereRadii{end+1}   = sphereRadii{1};
    sphereCenters{end+1} = sphereCenters{1};
    definingEdges{end+1} = definingEdges{1};


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

    % Find the modpoints for each vertex
    midpoints = (vertices(2:end,:)+vertices(1:end-1,:))/2;
    midpoints = [midpoints; midpoints(1,:)]; % Repeat first row, to last, to match how point is similarly repeated


    for this_vertex = 1:NumUniqueVerticies

        radiiFromVertexToEdge  = sphereRadii{this_vertex};
        sphereCenterArray      = sphereCenters{this_vertex};
        edgesConstrainingRadii = definingEdges{this_vertex};

        nexttile;

        grid on
        grid minor
        hold on
        axis equal


        % Plot the polytope in black dots connected by lines
        plot(vertices(:,1),vertices(:,2),'b.-','Linewidth',2, 'MarkerSize',10);

        % Plot the "ghostlines"
        unit_tangent_vectors = unit_normal_vectors*[0 1; -1 0];
        for ith_vertex = 1:length(vertices(:,1))-1
            ghostEnds = [...
                vertices(ith_vertex,:)+2*sizePlot*unit_tangent_vectors(ith_vertex,:);
                vertices(ith_vertex,:)-2*sizePlot*unit_tangent_vectors(ith_vertex,:);
                ];
            plot(ghostEnds(:,1),ghostEnds(:,2),'-','LineWidth',0.5,'Color',0.7*[1 1 1]);

        end

        % Plot the current vertex "ghostline"
        ghostEnds = [...
            vertices(this_vertex,:)+2*sizePlot*vertex_projection_vectors(this_vertex,:);
            vertices(this_vertex,:)-2*sizePlot*vertex_projection_vectors(this_vertex,:);
            ];
        plot(ghostEnds(:,1),ghostEnds(:,2),'-','LineWidth',0.5,'Color',0.7*[0 1 0]);


        % Label the vertices with their numbers
        for ith_vertex = 1:length(vertices(:,1))-1
            text(vertices(ith_vertex,1)+nudge, vertices(ith_vertex,2),...
                sprintf('%.0d',ith_vertex));
        end

        % Label the edges with their numbers
        for ith_edge = 1:length(vertices(:,1))-1
            text(midpoints(ith_edge,1)+nudge, midpoints(ith_edge,2),...
                sprintf('%.0d',ith_edge),'Color',[0 1 0]);
        end

        % Plot the spheres, and label their centers
        for ith_sphere = 1:length(radiiFromVertexToEdge)
            circleCenter = sphereCenterArray(ith_sphere,:);
            circleRadius = radiiFromVertexToEdge(ith_sphere,1);
            circleEdge = edgesConstrainingRadii(ith_sphere);

            % Plot the circle center as a large dot, and store the color
            h_fig = plot(circleCenter(1,1),circleCenter(1,2),'.','MarkerSize',20);
            colorUsed = get(h_fig,'Color');

            % Plot the circle boundary in same color
            fcn_geometry_plotCircle(circleCenter,circleRadius,colorUsed,fig_num);

            text(circleCenter(1,1)+nudge, circleCenter(1,2),...
                sprintf('%.0dto %.0d',this_vertex,circleEdge), 'Color',colorUsed);
        end

        % Make axis slightly larger?
        if flag_rescale_axis
            axis(goodAxis);
        end
        title(sprintf('Vertex: %.0d',this_vertex));
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

%% fcn_INTERNAL_isPointOnEdge
function isOnEdge = fcn_INTERNAL_isPointOnEdge(edgeStart,edgeEnd,testPoint)
methodToCheck = 1;
tolerance = 1E-5;

edgeLength     = real(sum((edgeEnd - edgeStart).^2,2).^0.5);

% For debugging

if 1==0
    figure(48484);
    clf
    inputdata = [edgeStart; edgeEnd];
    plot(inputdata(:,1),inputdata(:,2),'.-','MarkerSize',10,'LineWidth',3);
    hold on;
    axis equal
    plot(testPoint(:,1),testPoint(:,2),'.','MarkerSize',20);
end

if methodToCheck==1
    % Use vector projection
    edgeVector = edgeEnd - edgeStart;
    testVector = testPoint - edgeStart;
    unit_edgeVector = edgeVector./edgeLength;


    % Do dot product to see if within length
    dotProduct = sum(unit_edgeVector.*testVector,2);
    if dotProduct<(0-tolerance) || dotProduct>(edgeLength+tolerance)
        isOnEdge = 0;
        return;
    end

    % Do dot product to see if not orthogonal
    ortho_edgeVector = unit_edgeVector*[0 1; -1 0];
    dotProduct = sum(ortho_edgeVector.*testVector,2);
    if dotProduct<(0-tolerance) || dotProduct>(0+tolerance)
        isOnEdge = 0;
        return;
    end
    isOnEdge = 1;

else

    % Use the triangle inequality to check if test point is on edge
    pointDistances = sum((testPoint - edgeStart).^2,2).^0.5 + sum((testPoint - edgeEnd).^2,2).^0.5;
    offsetError = abs(edgeLength - pointDistances);
    if offsetError>tolerance
        isOnEdge = 0;
    else
        isOnEdge = 1;
    end
end

end % Ends fcn_INTERNAL_isPointOnEdge

%% fcn_INTERNAL_isPointWithinVertexProjections
function isWithin = fcn_INTERNAL_isPointWithinVertexProjections(...
    testPoint, vertexPoint, vertexProjection, ...
    currentVertexEdgeNormal, previousVertexEdgeNormal)


testVector = testPoint - vertexPoint;

currentEdgeProjection = currentVertexEdgeNormal*[0 1; -1 0];
previousEdgeProjection = previousVertexEdgeNormal*[0 1; -1 0];

% Make sure they are oriented correctly
if sum(currentEdgeProjection.*vertexProjection,2)<0
    currentEdgeProjection = currentEdgeProjection*(-1);
end
if sum(previousEdgeProjection.*vertexProjection,2)<0
    previousEdgeProjection = previousEdgeProjection*(-1);
end


% For debugging
if 1==0
    figure(77777);
    clf
    
    plot(vertexPoint(:,1),vertexPoint(:,2),'.-','MarkerSize',10,'LineWidth',3,'Color',[0 0 0]);
    hold on;
    grid on;
    axis equal
    plot(testPoint(:,1),testPoint(:,2),'.','MarkerSize',20,'Color',[0 1 0]);
    quiver(vertexPoint(1,1),vertexPoint(1,2),testVector(1,1),testVector(1,2),0,'Color',[0 0 0]);
    quiver(vertexPoint(1,1),vertexPoint(1,2),currentEdgeProjection(1,1),currentEdgeProjection(1,2),0,'Color',[0 1 0]);
    quiver(vertexPoint(1,1),vertexPoint(1,2),previousEdgeProjection(1,1),previousEdgeProjection(1,2),0,'Color',[1 0 0]);
    quiver(vertexPoint(1,1),vertexPoint(1,2),currentVertexEdgeNormal(1,1),currentVertexEdgeNormal(1,2),0,'Color',0.5*[0 1 0]);
    quiver(vertexPoint(1,1),vertexPoint(1,2),previousVertexEdgeNormal(1,1),previousVertexEdgeNormal(1,2),0,'Color',0.5*[1 0 0]);




end



dotProductCurrent = sum(currentVertexEdgeNormal.*testVector,2);
dotProductPrevious = sum(previousVertexEdgeNormal.*testVector,2);
projectionCurrent = sum(currentEdgeProjection.*testVector,2);
projectionPrevious = sum(previousEdgeProjection.*testVector,2);
conditionsToTest = [dotProductCurrent; dotProductPrevious; projectionCurrent; projectionPrevious];

isWithin = 1;
tolerance = 1E-6;
if any(conditionsToTest<(0-tolerance))
    isWithin = 0;
end
end % Ends fcn_INTERNAL_isPointWithinVertexProjections