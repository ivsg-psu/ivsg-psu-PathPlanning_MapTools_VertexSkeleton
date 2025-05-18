% script_test_fcn_VSkel_polytopeFindVertexSkeleton
% Tests function: fcn_VSkel_polytopeFindVertexSkeleton

% REVISION HISTORY:
% 2022_02_15
% -- first written by S. Brennan
% 2025_04_29 by Sean Brennan
% -- updated main script for VSkel library implementation

close all;


%% Demonstration Examples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  _____                                 _             _   _               ______                           _
% |  __ \                               | |           | | (_)             |  ____|                         | |
% | |  | | ___ _ __ ___   ___  _ __  ___| |_ _ __ __ _| |_ _  ___  _ __   | |__  __  ____ _ _ __ ___  _ __ | | ___  ___
% | |  | |/ _ \ '_ ` _ \ / _ \| '_ \/ __| __| '__/ _` | __| |/ _ \| '_ \  |  __| \ \/ / _` | '_ ` _ \| '_ \| |/ _ \/ __|
% | |__| |  __/ | | | | | (_) | | | \__ \ |_| | | (_| | |_| | (_) | | | | | |____ >  < (_| | | | | | | |_) | |  __/\__ \
% |_____/ \___|_| |_| |_|\___/|_| |_|___/\__|_|  \__,_|\__|_|\___/|_| |_| |______/_/\_\__,_|_| |_| |_| .__/|_|\___||___/
%                                                                                                    | |
%                                                                                                    |_|
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Demonstration%20Examples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Demonstration case 1: non-normal wall shrinking
fig_num = 9001;
figure(fig_num);
clf;

vertices = [0 0; 10 0; 5 10; 0 5; 0 0];
[cut_distance, vertexSkeleton] = fcn_VSkel_polytopeFindVertexSkeleton(vertices,fig_num);

% Check variable types
assert(isnumeric(cut_distance));
assert(isstruct(vertexSkeleton));

% Check variable sizes
assert(length(vertexSkeleton)==length(cut_distance));

for ith_depth = 1:length(cut_distance)    
    for ith_polytope = 1:length(vertexSkeleton(ith_depth).polytope)
        Nvertices = length(vertexSkeleton(ith_depth).polytope(ith_polytope).vertices(:,1));
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).vector_direction_of_unit_cut(:,1))==Nvertices);
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).unit_vertex_projection_vectors(:,1))==Nvertices);
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).flag_vertexIsNonConvex(:,1))==Nvertices);
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).intersection_points(:,1))==Nvertices);
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).min_cut(:,1))==1);
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).indices_repeated(:,1))==length(vertexSkeleton(ith_depth).polytope(ith_polytope).boundaryEngagedAtMinCut(:,1)));

    end
end

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

% Make a plot of nested cuts?
if 1==1
    figure(4747);
    grid on
    grid minor
    hold on
    axis equal

    for cut = 0:0.5:cut_distance(end)

        % Find the shape that is less than or equal to the cut
        shape_index = find(cut_distance<=cut,1,'last');

        % Grab vertices to start from, cut to start from
        starting_vertices = vertexSkeleton(shape_index).polytope(1).vertices;
        new_projection_vectors = vertexSkeleton(shape_index).polytope(1).vector_direction_of_unit_cut;
        starting_cut      = cut_distance(shape_index);

        % Calculate projection distance
        projection_distance = cut - starting_cut;

        % Determine final vertices
        final_vertices = starting_vertices + new_projection_vectors*projection_distance;

        % Plot results
        plot(final_vertices(:,1),final_vertices(:,2),'.-','Linewidth',2,'Markersize',20);
    end
end

%% Basic testing examples - Line Segments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  ____            _        _______        _   _                ______                           _
% |  _ \          (_)      |__   __|      | | (_)              |  ____|                         | |
% | |_) | __ _ ___ _  ___     | | ___  ___| |_ _ _ __   __ _   | |__  __  ____ _ _ __ ___  _ __ | | ___  ___
% |  _ < / _` / __| |/ __|    | |/ _ \/ __| __| | '_ \ / _` |  |  __| \ \/ / _` | '_ ` _ \| '_ \| |/ _ \/ __|
% | |_) | (_| \__ \ | (__     | |  __/\__ \ |_| | | | | (_| |  | |____ >  < (_| | | | | | | |_) | |  __/\__ \
% |____/ \__,_|___/_|\___|    |_|\___||___/\__|_|_| |_|\__, |  |______/_/\_\__,_|_| |_| |_| .__/|_|\___||___/
%                                                       __/ |                             | |
%                                                      |___/                              |_|
%      __     _      _               _____                                 _
%      \ \   | |    (_)             / ____|                               | |
%  _____\ \  | |     _ _ __   ___  | (___   ___  __ _ _ __ ___   ___ _ __ | |_ ___
% |______> > | |    | | '_ \ / _ \  \___ \ / _ \/ _` | '_ ` _ \ / _ \ '_ \| __/ __|
%       / /  | |____| | | | |  __/  ____) |  __/ (_| | | | | | |  __/ | | | |_\__ \
%      /_/   |______|_|_| |_|\___| |_____/ \___|\__, |_| |_| |_|\___|_| |_|\__|___/
%                                                __/ |
%                                               |___/
% See: http://patorjk.com/software/taag/#p=display&v=0&f=Big&t=Basic%20Testing%20%20Examples%20%0A-%3E%20Line%20Segments%0A
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Basic example of vertex calculation - three point test, horizontal line segment
fig_num = 1001;
figure(fig_num);
clf;

vertices = [0 0; 2 0; 0 0];
[cut_distance, vertexSkeleton] = fcn_VSkel_polytopeFindVertexSkeleton(vertices,fig_num);

% Check variable types
assert(isnumeric(cut_distance));
assert(isstruct(vertexSkeleton));

% Check variable sizes
assert(length(vertexSkeleton)==length(cut_distance));

for ith_depth = 1:length(cut_distance)    
    for ith_polytope = 1:length(vertexSkeleton(ith_depth).polytope)
        Nvertices = length(vertexSkeleton(ith_depth).polytope(ith_polytope).vertices(:,1));
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).vector_direction_of_unit_cut(:,1))==Nvertices);
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).unit_vertex_projection_vectors(:,1))==Nvertices);
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).flag_vertexIsNonConvex(:,1))==Nvertices);
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).intersection_points(:,1))==Nvertices);
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).min_cut(:,1))==1);
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).indices_repeated(:,1))==length(vertexSkeleton(ith_depth).polytope(ith_polytope).boundaryEngagedAtMinCut(:,1)));

    end
end

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Basic example of vertex calculation - three point test, vertical line segment
fig_num = 1002;
figure(fig_num);
clf;

vertices = [0 0; 0 2; 0 0];
[cut_distance, vertexSkeleton] = fcn_VSkel_polytopeFindVertexSkeleton(vertices,fig_num);

% Check variable types
assert(isnumeric(cut_distance));
assert(isstruct(vertexSkeleton));

% Check variable sizes
assert(length(vertexSkeleton)==length(cut_distance));

for ith_depth = 1:length(cut_distance)    
    for ith_polytope = 1:length(vertexSkeleton(ith_depth).polytope)
        Nvertices = length(vertexSkeleton(ith_depth).polytope(ith_polytope).vertices(:,1));
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).vector_direction_of_unit_cut(:,1))==Nvertices);
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).unit_vertex_projection_vectors(:,1))==Nvertices);
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).flag_vertexIsNonConvex(:,1))==Nvertices);
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).intersection_points(:,1))==Nvertices);
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).min_cut(:,1))==1);
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).indices_repeated(:,1))==length(vertexSkeleton(ith_depth).polytope(ith_polytope).boundaryEngagedAtMinCut(:,1)));

    end
end

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Basic example of vertex calculation - three point test, arbitrary line segment
fig_num = 1003;
figure(fig_num);
clf;

vertices = [-1 2; 3 4; -1 2];
[cut_distance, vertexSkeleton] = fcn_VSkel_polytopeFindVertexSkeleton(vertices,fig_num);

% Check variable types
assert(isnumeric(cut_distance));
assert(isstruct(vertexSkeleton));

% Check variable sizes
assert(length(vertexSkeleton)==length(cut_distance));

for ith_depth = 1:length(cut_distance)    
    for ith_polytope = 1:length(vertexSkeleton(ith_depth).polytope)
        Nvertices = length(vertexSkeleton(ith_depth).polytope(ith_polytope).vertices(:,1));
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).vector_direction_of_unit_cut(:,1))==Nvertices);
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).unit_vertex_projection_vectors(:,1))==Nvertices);
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).flag_vertexIsNonConvex(:,1))==Nvertices);
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).intersection_points(:,1))==Nvertices);
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).min_cut(:,1))==1);
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).indices_repeated(:,1))==length(vertexSkeleton(ith_depth).polytope(ith_polytope).boundaryEngagedAtMinCut(:,1)));

    end
end

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));



%% Basic testing examples - Symmetric 2D Polytopes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  ____            _        _______        _   _                ______                           _
% |  _ \          (_)      |__   __|      | | (_)              |  ____|                         | |
% | |_) | __ _ ___ _  ___     | | ___  ___| |_ _ _ __   __ _   | |__  __  ____ _ _ __ ___  _ __ | | ___  ___
% |  _ < / _` / __| |/ __|    | |/ _ \/ __| __| | '_ \ / _` |  |  __| \ \/ / _` | '_ ` _ \| '_ \| |/ _ \/ __|
% | |_) | (_| \__ \ | (__     | |  __/\__ \ |_| | | | | (_| |  | |____ >  < (_| | | | | | | |_) | |  __/\__ \
% |____/ \__,_|___/_|\___|    |_|\___||___/\__|_|_| |_|\__, |  |______/_/\_\__,_|_| |_| |_| .__/|_|\___||___/
%                                                       __/ |                             | |
%                                                      |___/                              |_|
%      __      _____                                _        _        ___  _____    _____      _       _
%      \ \    / ____|                              | |      (_)      |__ \|  __ \  |  __ \    | |     | |
%  _____\ \  | (___  _   _ _ __ ___  _ __ ___   ___| |_ _ __ _  ___     ) | |  | | | |__) |__ | |_   _| |_ ___  _ __   ___  ___
% |______> >  \___ \| | | | '_ ` _ \| '_ ` _ \ / _ \ __| '__| |/ __|   / /| |  | | |  ___/ _ \| | | | | __/ _ \| '_ \ / _ \/ __|
%       / /   ____) | |_| | | | | | | | | | | |  __/ |_| |  | | (__   / /_| |__| | | |  | (_) | | |_| | || (_) | |_) |  __/\__ \
%      /_/   |_____/ \__, |_| |_| |_|_| |_| |_|\___|\__|_|  |_|\___| |____|_____/  |_|   \___/|_|\__, |\__\___/| .__/ \___||___/
%                     __/ |                                                                       __/ |        | |
%                    |___/                                                                       |___/         |_|
% See: http://patorjk.com/software/taag/#p=display&v=0&f=Big&t=Basic%20Testing%20%20Examples%20%0A-%3E%20Symmetric%202D%20Polytopes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Example case 1: symmetric tall triangle
fig_num = 2001;
figure(fig_num);
clf;

vertices = [0 0; 10 0; 5 10; 0 0];
[cut_distance, vertexSkeleton] = fcn_VSkel_polytopeFindVertexSkeleton(vertices,fig_num);

% Check variable types
assert(isnumeric(cut_distance));
assert(isstruct(vertexSkeleton));

% Check variable sizes
assert(length(vertexSkeleton)==length(cut_distance));

for ith_depth = 1:length(cut_distance)    
    for ith_polytope = 1:length(vertexSkeleton(ith_depth).polytope)
        Nvertices = length(vertexSkeleton(ith_depth).polytope(ith_polytope).vertices(:,1));
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).vector_direction_of_unit_cut(:,1))==Nvertices);
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).unit_vertex_projection_vectors(:,1))==Nvertices);
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).flag_vertexIsNonConvex(:,1))==Nvertices);
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).intersection_points(:,1))==Nvertices);
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).min_cut(:,1))==1);
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).indices_repeated(:,1))==length(vertexSkeleton(ith_depth).polytope(ith_polytope).boundaryEngagedAtMinCut(:,1)));

    end
end

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Example case 2: square - will have 2 nested solutions
fig_num = 2002;
figure(fig_num);
clf;

vertices = [0 0; 10 0; 10 10; 0 10; 0 0];
[cut_distance, vertexSkeleton] = fcn_VSkel_polytopeFindVertexSkeleton(vertices,fig_num);

% Check variable types
assert(isnumeric(cut_distance));
assert(isstruct(vertexSkeleton));

% Check variable sizes
assert(length(vertexSkeleton)==length(cut_distance));

for ith_depth = 1:length(cut_distance)    
    for ith_polytope = 1:length(vertexSkeleton(ith_depth).polytope)
        Nvertices = length(vertexSkeleton(ith_depth).polytope(ith_polytope).vertices(:,1));
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).vector_direction_of_unit_cut(:,1))==Nvertices);
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).unit_vertex_projection_vectors(:,1))==Nvertices);
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).flag_vertexIsNonConvex(:,1))==Nvertices);
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).intersection_points(:,1))==Nvertices);
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).min_cut(:,1))==1);
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).indices_repeated(:,1))==length(vertexSkeleton(ith_depth).polytope(ith_polytope).boundaryEngagedAtMinCut(:,1)));

    end
end

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Example case 3: wide rectangle
fig_num = 2003;
figure(fig_num);
clf;

vertices = [0 0; 10 0; 10 5; 0 5; 0 0];
[cut_distance, vertexSkeleton] = fcn_VSkel_polytopeFindVertexSkeleton(vertices,fig_num);

% Check variable types
assert(isnumeric(cut_distance));
assert(isstruct(vertexSkeleton));

% Check variable sizes
assert(length(vertexSkeleton)==length(cut_distance));

for ith_depth = 1:length(cut_distance)    
    for ith_polytope = 1:length(vertexSkeleton(ith_depth).polytope)
        Nvertices = length(vertexSkeleton(ith_depth).polytope(ith_polytope).vertices(:,1));
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).vector_direction_of_unit_cut(:,1))==Nvertices);
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).unit_vertex_projection_vectors(:,1))==Nvertices);
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).flag_vertexIsNonConvex(:,1))==Nvertices);
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).intersection_points(:,1))==Nvertices);
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).min_cut(:,1))==1);
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).indices_repeated(:,1))==length(vertexSkeleton(ith_depth).polytope(ith_polytope).boundaryEngagedAtMinCut(:,1)));

    end
end

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));
%% Example case 4: tall rectangle
fig_num = 2004;
figure(fig_num);
clf;

vertices = [0 0; 0.5 0; 0.5 1; 0 1; 0 0]*10;
[cut_distance, vertexSkeleton] = fcn_VSkel_polytopeFindVertexSkeleton(vertices,fig_num);

% Check variable types
assert(isnumeric(cut_distance));
assert(isstruct(vertexSkeleton));

% Check variable sizes
assert(length(vertexSkeleton)==length(cut_distance));

for ith_depth = 1:length(cut_distance)    
    for ith_polytope = 1:length(vertexSkeleton(ith_depth).polytope)
        Nvertices = length(vertexSkeleton(ith_depth).polytope(ith_polytope).vertices(:,1));
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).vector_direction_of_unit_cut(:,1))==Nvertices);
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).unit_vertex_projection_vectors(:,1))==Nvertices);
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).flag_vertexIsNonConvex(:,1))==Nvertices);
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).intersection_points(:,1))==Nvertices);
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).min_cut(:,1))==1);
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).indices_repeated(:,1))==length(vertexSkeleton(ith_depth).polytope(ith_polytope).boundaryEngagedAtMinCut(:,1)));

    end
end

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Basic testing examples - Assymetric 2D Convex Polytopes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  ____            _        _______        _   _                ______                           _
% |  _ \          (_)      |__   __|      | | (_)              |  ____|                         | |
% | |_) | __ _ ___ _  ___     | | ___  ___| |_ _ _ __   __ _   | |__  __  ____ _ _ __ ___  _ __ | | ___  ___
% |  _ < / _` / __| |/ __|    | |/ _ \/ __| __| | '_ \ / _` |  |  __| \ \/ / _` | '_ ` _ \| '_ \| |/ _ \/ __|
% | |_) | (_| \__ \ | (__     | |  __/\__ \ |_| | | | | (_| |  | |____ >  < (_| | | | | | | |_) | |  __/\__ \
% |____/ \__,_|___/_|\___|    |_|\___||___/\__|_|_| |_|\__, |  |______/_/\_\__,_|_| |_| |_| .__/|_|\___||___/
%                                                       __/ |                             | |
%                                                      |___/                              |_|
%      __                                           _        _        ___  _____     _____                            _____      _       _
%      \ \       /\                                | |      (_)      |__ \|  __ \   / ____|                          |  __ \    | |     | |
%  _____\ \     /  \   ___ ___ _   _ _ __ ___   ___| |_ _ __ _  ___     ) | |  | | | |     ___  _ ____   _______  __ | |__) |__ | |_   _| |_ ___  _ __   ___  ___
% |______> >   / /\ \ / __/ __| | | | '_ ` _ \ / _ \ __| '__| |/ __|   / /| |  | | | |    / _ \| '_ \ \ / / _ \ \/ / |  ___/ _ \| | | | | __/ _ \| '_ \ / _ \/ __|
%       / /   / ____ \\__ \__ \ |_| | | | | | |  __/ |_| |  | | (__   / /_| |__| | | |___| (_) | | | \ V /  __/>  <  | |  | (_) | | |_| | || (_) | |_) |  __/\__ \
%      /_/   /_/    \_\___/___/\__, |_| |_| |_|\___|\__|_|  |_|\___| |____|_____/   \_____\___/|_| |_|\_/ \___/_/\_\ |_|   \___/|_|\__, |\__\___/| .__/ \___||___/
%                               __/ |                                                                                               __/ |        | |
%                              |___/                                                                                               |___/         |_|
% See: http://patorjk.com/software/taag/#p=display&v=0&f=Big&t=Basic%20Testing%20%20Examples%20%0A-%3E%20Assymetric%202D%20Convex%20Polytopes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Example case 1: goofy polytope
fig_num = 3001;
figure(fig_num);
clf;

vertices = [0 0; 10 0; 5 15; 4 17; 1 13; 0 5; 0 0];
[cut_distance, vertexSkeleton] = fcn_VSkel_polytopeFindVertexSkeleton(vertices,fig_num);

% Check variable types
assert(isnumeric(cut_distance));
assert(isstruct(vertexSkeleton));

% Check variable sizes
assert(length(vertexSkeleton)==length(cut_distance));

for ith_depth = 1:length(cut_distance)    
    for ith_polytope = 1:length(vertexSkeleton(ith_depth).polytope)
        Nvertices = length(vertexSkeleton(ith_depth).polytope(ith_polytope).vertices(:,1));
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).vector_direction_of_unit_cut(:,1))==Nvertices);
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).unit_vertex_projection_vectors(:,1))==Nvertices);
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).flag_vertexIsNonConvex(:,1))==Nvertices);
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).intersection_points(:,1))==Nvertices);
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).min_cut(:,1))==1);
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).indices_repeated(:,1))==length(vertexSkeleton(ith_depth).polytope(ith_polytope).boundaryEngagedAtMinCut(:,1)));

    end
end

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));


%% Random polytope calculation
fig_num = 3002;
figure(fig_num);
clf;

rng(1); % Set the random number generator

% Set up polytopes
load('testData_fcn_VSkel_plotPolytopes.mat','polytopes','polytopes2');

% Pick a random polytope
Npolys = length(polytopes);
rand_poly = 1+floor(rand*Npolys);
shrinker = polytopes(rand_poly);
vertices = shrinker.vertices;

[cut_distance, vertexSkeleton] = fcn_VSkel_polytopeFindVertexSkeleton(vertices,fig_num);

% Check variable types
assert(isnumeric(cut_distance));
assert(isstruct(vertexSkeleton));

% Check variable sizes
assert(length(vertexSkeleton)==length(cut_distance));

for ith_depth = 1:length(cut_distance)    
    for ith_polytope = 1:length(vertexSkeleton(ith_depth).polytope)
        Nvertices = length(vertexSkeleton(ith_depth).polytope(ith_polytope).vertices(:,1));
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).vector_direction_of_unit_cut(:,1))==Nvertices);
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).unit_vertex_projection_vectors(:,1))==Nvertices);
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).flag_vertexIsNonConvex(:,1))==Nvertices);
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).intersection_points(:,1))==Nvertices);
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).min_cut(:,1))==1);
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).indices_repeated(:,1))==length(vertexSkeleton(ith_depth).polytope(ith_polytope).boundaryEngagedAtMinCut(:,1)));

    end
end

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

% Make a plot of nested cuts?
if 1==1
    figure(4747);
    clf;
    C = colororder; 
    Ncolors = length(C(:,1));

    grid on
    grid minor
    hold on
    axis equal

    edge_cut_step = 0.002;
    for cut = edge_cut_step:edge_cut_step:(shrinker.max_radius+edge_cut_step)

        % Find the shape that is less than or equal to the cut
        shape_index = find(cut_distance<=cut,1,'last');

        % Grab vertices to start from, cut to start from
        starting_vertices = vertexSkeleton(shape_index).polytope(1).vertices;
        new_projection_vectors = vertexSkeleton(shape_index).polytope(1).vector_direction_of_unit_cut;
        starting_cut      = cut_distance(shape_index);

        % Calculate projection distance
        projection_distance = cut - starting_cut;

        % Determine final vertices
        final_vertices = starting_vertices + new_projection_vectors*projection_distance;

        % Find the plotting color
 
        current_cut_color_index = mod(shape_index-1,Ncolors)+1;
        current_cut_base_color = C(current_cut_color_index,:);
        if 1 == 0
            % Blend colors?
            next_color_index = mod(shape_index,Ncolors)+1;
            next_cut_base_color    = C(next_color_index,:);

            next_cut_index    = min(shape_index+1,length(cut_distance(:,1)));
            next_cut = cut_distance(next_cut_index);

            percent_mix = projection_distance/(next_cut-starting_cut);
            if percent_mix>1
                percent_mix = 1;
            end

            current_color = current_cut_base_color*(1-percent_mix) + next_cut_base_color*(percent_mix);
        else
            current_color = current_cut_base_color;
        end

        % Plot results
        plot(final_vertices(:,1),final_vertices(:,2),'.-','Linewidth',2,'Markersize',20, 'Color',current_color);
    end
end

%% Basic testing examples - Assymetric 2D Non-Convex Polytopes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  ____            _        _______        _   _                ______                           _
% |  _ \          (_)      |__   __|      | | (_)              |  ____|                         | |
% | |_) | __ _ ___ _  ___     | | ___  ___| |_ _ _ __   __ _   | |__  __  ____ _ _ __ ___  _ __ | | ___  ___
% |  _ < / _` / __| |/ __|    | |/ _ \/ __| __| | '_ \ / _` |  |  __| \ \/ / _` | '_ ` _ \| '_ \| |/ _ \/ __|
% | |_) | (_| \__ \ | (__     | |  __/\__ \ |_| | | | | (_| |  | |____ >  < (_| | | | | | | |_) | |  __/\__ \
% |____/ \__,_|___/_|\___|    |_|\___||___/\__|_|_| |_|\__, |  |______/_/\_\__,_|_| |_| |_| .__/|_|\___||___/
%                                                       __/ |                             | |
%                                                      |___/                              |_|
%      __                                           _        _        ___  _____    _   _                    _____                            _____      _       _
%      \ \       /\                                | |      (_)      |__ \|  __ \  | \ | |                  / ____|                          |  __ \    | |     | |
%  _____\ \     /  \   ___ ___ _   _ _ __ ___   ___| |_ _ __ _  ___     ) | |  | | |  \| | ___  _ __ ______| |     ___  _ ____   _______  __ | |__) |__ | |_   _| |_ ___  _ __   ___  ___
% |______> >   / /\ \ / __/ __| | | | '_ ` _ \ / _ \ __| '__| |/ __|   / /| |  | | | . ` |/ _ \| '_ \______| |    / _ \| '_ \ \ / / _ \ \/ / |  ___/ _ \| | | | | __/ _ \| '_ \ / _ \/ __|
%       / /   / ____ \\__ \__ \ |_| | | | | | |  __/ |_| |  | | (__   / /_| |__| | | |\  | (_) | | | |     | |___| (_) | | | \ V /  __/>  <  | |  | (_) | | |_| | || (_) | |_) |  __/\__ \
%      /_/   /_/    \_\___/___/\__, |_| |_| |_|\___|\__|_|  |_|\___| |____|_____/  |_| \_|\___/|_| |_|      \_____\___/|_| |_|\_/ \___/_/\_\ |_|   \___/|_|\__, |\__\___/| .__/ \___||___/
%                               __/ |                                                                                                                       __/ |        | |
%                              |___/                                                                                                                       |___/         |_|
% See: http://patorjk.com/software/taag/#p=display&v=0&f=Big&t=Basic%20Testing%20%20Examples%20%0A-%3E%20Assymetric%202D%20Non-Convex%20Polytopes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Basic example of vertex calculation - non-convex 2D polytope
fig_num = 4001;
figure(fig_num);
clf;

% this polytope has a vertical wall
vertices = [0 0; 3 0; 5 5; 7 0; 10 0; 5 10; 0 5; 0 0];
[cut_distance, vertexSkeleton] = fcn_VSkel_polytopeFindVertexSkeleton(vertices,fig_num);

% Check variable types
assert(isnumeric(cut_distance));
assert(isstruct(vertexSkeleton));

% Check variable sizes
assert(length(vertexSkeleton)==length(cut_distance));

for ith_depth = 1:length(cut_distance)    
    for ith_polytope = 1:length(vertexSkeleton(ith_depth).polytope)
        Nvertices = length(vertexSkeleton(ith_depth).polytope(ith_polytope).vertices(:,1));
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).vector_direction_of_unit_cut(:,1))==Nvertices);
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).unit_vertex_projection_vectors(:,1))==Nvertices);
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).flag_vertexIsNonConvex(:,1))==Nvertices);
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).intersection_points(:,1))==Nvertices);
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).min_cut(:,1))==1);
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).indices_repeated(:,1))==length(vertexSkeleton(ith_depth).polytope(ith_polytope).boundaryEngagedAtMinCut(:,1)));

    end
end

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

% Make a plot of nested cuts?
if 1==1
    figure(4747);
    clf;
    C = colororder; 
    Ncolors = length(C(:,1));

    grid on
    grid minor
    hold on
    axis equal

    edge_cut_step = 0.05;
    for cut = 0:edge_cut_step:cut_distance(end)

        % Find the shape that is less than or equal to the cut
        shape_index = find(cut_distance<=cut,1,'last');

        Npolys = length(vertexSkeleton(shape_index).polytope);
        for ith_polytope = 1:Npolys
            % Grab vertices to start from, cut to start from
            starting_vertices = vertexSkeleton(shape_index).polytope(ith_polytope).vertices;
            new_projection_vectors = vertexSkeleton(shape_index).polytope(ith_polytope).vector_direction_of_unit_cut;
            starting_cut      = cut_distance(shape_index);

            % Calculate projection distance
            projection_distance = cut - starting_cut;

            % Determine final vertices
            final_vertices = starting_vertices + new_projection_vectors*projection_distance;

            % Find the plotting color

            current_cut_color_index = mod(shape_index-1,Ncolors)+1;
            current_cut_base_color = C(current_cut_color_index,:);
            if 1 == 0
                % Blend colors?
                next_color_index = mod(shape_index,Ncolors)+1;
                next_cut_base_color    = C(next_color_index,:);

                next_cut_index    = min(shape_index+1,length(cut_distance(:,1)));
                next_cut = cut_distance(next_cut_index);

                percent_mix = projection_distance/(next_cut-starting_cut);
                if percent_mix>1
                    percent_mix = 1;
                end

                current_color = current_cut_base_color*(1-percent_mix) + next_cut_base_color*(percent_mix);
            else
                current_color = current_cut_base_color;
            end

            % Plot results
            plot(final_vertices(:,1),final_vertices(:,2),'-','Linewidth',2,'Markersize',20, 'Color',current_color);
        end
    end
end


%% Basic example of vertex calculation - hard non-convex 2D polytope
fig_num =4002;
figure(fig_num);
clf;

% this polytope has a vertical wall
vertices = [0 0; 5 0; 6 4; 7 0; 10 0; 10 10; 4 10; 4 5; 0 5; 0 0];
[cut_distance, vertexSkeleton] = fcn_VSkel_polytopeFindVertexSkeleton(vertices,fig_num);

% Check variable types
assert(isnumeric(cut_distance));
assert(isstruct(vertexSkeleton));

% Check variable sizes
assert(length(vertexSkeleton)==length(cut_distance));

for ith_depth = 1:length(cut_distance)    
    for ith_polytope = 1:length(vertexSkeleton(ith_depth).polytope)
        Nvertices = length(vertexSkeleton(ith_depth).polytope(ith_polytope).vertices(:,1));
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).vector_direction_of_unit_cut(:,1))==Nvertices);
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).unit_vertex_projection_vectors(:,1))==Nvertices);
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).flag_vertexIsNonConvex(:,1))==Nvertices);
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).intersection_points(:,1))==Nvertices);
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).min_cut(:,1))==1);
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).indices_repeated(:,1))==length(vertexSkeleton(ith_depth).polytope(ith_polytope).boundaryEngagedAtMinCut(:,1)));

    end
end

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));


%% Basic example of vertex calculation - symmetric non-convex vertices
fig_num = 4003;
figure(fig_num);
clf;

% this polytope has two nonconvex vertices facing each other
vertices = [5 4; 6 0; 10 0; 10 10; 6 10; 5 6; 4 10; 0 10; 0 0; 4 0; 5 4];
[cut_distance, vertexSkeleton] = fcn_VSkel_polytopeFindVertexSkeleton(vertices,fig_num);

% Check variable types
assert(isnumeric(cut_distance));
assert(isstruct(vertexSkeleton));

% Check variable sizes
assert(length(vertexSkeleton)==length(cut_distance));

for ith_depth = 1:length(cut_distance)    
    for ith_polytope = 1:length(vertexSkeleton(ith_depth).polytope)
        Nvertices = length(vertexSkeleton(ith_depth).polytope(ith_polytope).vertices(:,1));
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).vector_direction_of_unit_cut(:,1))==Nvertices);
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).unit_vertex_projection_vectors(:,1))==Nvertices);
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).flag_vertexIsNonConvex(:,1))==Nvertices);
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).intersection_points(:,1))==Nvertices);
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).min_cut(:,1))==1);
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).indices_repeated(:,1))==length(vertexSkeleton(ith_depth).polytope(ith_polytope).boundaryEngagedAtMinCut(:,1)));

    end
end

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Basic example of vertex calculation - non-symmetric non-convex vertices
fig_num = 4004;
figure(fig_num);
clf;

% this polytope has two nonconvex vertices facing each other
vertices = [5 4; 6 0; 10 0; 10 7; 8 8; 10 9; 10 15; 0 15; 0 0; 4 0; 5 4];
[cut_distance, vertexSkeleton] = fcn_VSkel_polytopeFindVertexSkeleton(vertices,fig_num);

% Check variable types
assert(isnumeric(cut_distance));
assert(isstruct(vertexSkeleton));

% Check variable sizes
assert(length(vertexSkeleton)==length(cut_distance));

for ith_depth = 1:length(cut_distance)    
    for ith_polytope = 1:length(vertexSkeleton(ith_depth).polytope)
        Nvertices = length(vertexSkeleton(ith_depth).polytope(ith_polytope).vertices(:,1));
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).vector_direction_of_unit_cut(:,1))==Nvertices);
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).unit_vertex_projection_vectors(:,1))==Nvertices);
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).flag_vertexIsNonConvex(:,1))==Nvertices);
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).intersection_points(:,1))==Nvertices);
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).min_cut(:,1))==1);
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).indices_repeated(:,1))==length(vertexSkeleton(ith_depth).polytope(ith_polytope).boundaryEngagedAtMinCut(:,1)));

    end
end

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Basic example of vertex calculation - non-symmetric non-convex vertices
fig_num = 4005;
figure(fig_num);
clf;

% this polytope has two nonconvex vertices facing each other
vertices = [5 4; 6 0; 10 0; 10 7; 8 8; 10 9; 15 15; 0 15; 0 0; 4 0; 5 4];
[cut_distance, vertexSkeleton] = fcn_VSkel_polytopeFindVertexSkeleton(vertices,fig_num);

% Check variable types
assert(isnumeric(cut_distance));
assert(isstruct(vertexSkeleton));

% Check variable sizes
assert(length(vertexSkeleton)==length(cut_distance));

for ith_depth = 1:length(cut_distance)    
    for ith_polytope = 1:length(vertexSkeleton(ith_depth).polytope)
        Nvertices = length(vertexSkeleton(ith_depth).polytope(ith_polytope).vertices(:,1));
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).vector_direction_of_unit_cut(:,1))==Nvertices);
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).unit_vertex_projection_vectors(:,1))==Nvertices);
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).flag_vertexIsNonConvex(:,1))==Nvertices);
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).intersection_points(:,1))==Nvertices);
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).min_cut(:,1))==1);
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).indices_repeated(:,1))==length(vertexSkeleton(ith_depth).polytope(ith_polytope).boundaryEngagedAtMinCut(:,1)));

    end
end

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));


%% Basic testing examples - Enclosing 2D Polytopes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  ____            _        _______        _   _                ______                           _
% |  _ \          (_)      |__   __|      | | (_)              |  ____|                         | |
% | |_) | __ _ ___ _  ___     | | ___  ___| |_ _ _ __   __ _   | |__  __  ____ _ _ __ ___  _ __ | | ___  ___
% |  _ < / _` / __| |/ __|    | |/ _ \/ __| __| | '_ \ / _` |  |  __| \ \/ / _` | '_ ` _ \| '_ \| |/ _ \/ __|
% | |_) | (_| \__ \ | (__     | |  __/\__ \ |_| | | | | (_| |  | |____ >  < (_| | | | | | | |_) | |  __/\__ \
% |____/ \__,_|___/_|\___|    |_|\___||___/\__|_|_| |_|\__, |  |______/_/\_\__,_|_| |_| |_| .__/|_|\___||___/
%                                                       __/ |                             | |
%                                                      |___/                              |_|
%      __     ______            _           _               ___  _____     _____      _       _
%      \ \   |  ____|          | |         (_)             |__ \|  __ \   |  __ \    | |     | |
%  _____\ \  | |__   _ __   ___| | ___  ___ _ _ __   __ _     ) | |  | |  | |__) |__ | |_   _| |_ ___  _ __   ___  ___
% |______> > |  __| | '_ \ / __| |/ _ \/ __| | '_ \ / _` |   / /| |  | |  |  ___/ _ \| | | | | __/ _ \| '_ \ / _ \/ __|
%       / /  | |____| | | | (__| | (_) \__ \ | | | | (_| |  / /_| |__| |  | |  | (_) | | |_| | || (_) | |_) |  __/\__ \
%      /_/   |______|_| |_|\___|_|\___/|___/_|_| |_|\__, | |____|_____/   |_|   \___/|_|\__, |\__\___/| .__/ \___||___/
%                                                    __/ |                               __/ |        | |
%                                                   |___/                               |___/         |_|
% See: http://patorjk.com/software/taag/#p=display&v=0&f=Big&t=Basic%20Testing%20%20Examples%20%0A-%3E%20Enclosing%202D%20%20Polytopes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Example case 3: wide rectangle
fig_num = 5001;
figure(fig_num);
clf;

clear vertices
vertices{1} = [0 0; 10 0; 10 5; 0 5; 0 0];
vertices{2} = [2 2; 8 2; 8 3; 2 3; 2 2];
[cut_distance, vertexSkeleton] = fcn_VSkel_polytopeFindVertexSkeleton(vertices,fig_num);

% Check variable types
assert(isnumeric(cut_distance));
assert(isstruct(vertexSkeleton));

% Check variable sizes
assert(length(vertexSkeleton)==length(cut_distance));

for ith_depth = 1:length(cut_distance)    
    for ith_polytope = 1:length(vertexSkeleton(ith_depth).polytope)
        Nvertices = length(vertexSkeleton(ith_depth).polytope(ith_polytope).vertices(:,1));
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).vector_direction_of_unit_cut(:,1))==Nvertices);
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).unit_vertex_projection_vectors(:,1))==Nvertices);
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).flag_vertexIsNonConvex(:,1))==Nvertices);
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).intersection_points(:,1))==Nvertices);
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).min_cut(:,1))==1);
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).indices_repeated(:,1))==length(vertexSkeleton(ith_depth).polytope(ith_polytope).boundaryEngagedAtMinCut(:,1)));

    end
end

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Fast Mode Tests
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  ______        _     __  __           _        _______        _
% |  ____|      | |   |  \/  |         | |      |__   __|      | |
% | |__ __ _ ___| |_  | \  / | ___   __| | ___     | | ___  ___| |_ ___
% |  __/ _` / __| __| | |\/| |/ _ \ / _` |/ _ \    | |/ _ \/ __| __/ __|
% | | | (_| \__ \ |_  | |  | | (_) | (_| |  __/    | |  __/\__ \ |_\__ \
% |_|  \__,_|___/\__| |_|  |_|\___/ \__,_|\___|    |_|\___||___/\__|___/
%
%
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Fast%20Mode%20Tests
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Basic example of vertex calculation - non-normal wall shrinking, NO FIGURE
fig_num = 9901;
figure(fig_num);
close(fig_num);

vertices = [0 0; 10 0; 5 10; 0 5; 0 0];
[cut_distance, vertexSkeleton] = fcn_VSkel_polytopeFindVertexSkeleton(vertices,[]);

% Check variable types
assert(isnumeric(cut_distance));
assert(isstruct(vertexSkeleton));

% Check variable sizes
assert(length(vertexSkeleton)==length(cut_distance));

for ith_depth = 1:length(cut_distance)    
    for ith_polytope = 1:length(vertexSkeleton(ith_depth).polytope)
        Nvertices = length(vertexSkeleton(ith_depth).polytope(ith_polytope).vertices(:,1));
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).vector_direction_of_unit_cut(:,1))==Nvertices);
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).unit_vertex_projection_vectors(:,1))==Nvertices);
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).flag_vertexIsNonConvex(:,1))==Nvertices);
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).intersection_points(:,1))==Nvertices);
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).min_cut(:,1))==1);
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).indices_repeated(:,1))==length(vertexSkeleton(ith_depth).polytope(ith_polytope).boundaryEngagedAtMinCut(:,1)));

    end
end

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));

%% Basic example of vertex calculation - non-normal wall shrinking, NO FIGURE, FAST MODE
fig_num = 9902;
figure(fig_num);
close(fig_num);

vertices = [0 0; 10 0; 5 10; 0 5; 0 0];
[cut_distance, vertexSkeleton] = fcn_VSkel_polytopeFindVertexSkeleton(vertices,-1);

% Check variable types
assert(isnumeric(cut_distance));
assert(isstruct(vertexSkeleton));

% Check variable sizes
assert(length(vertexSkeleton)==length(cut_distance));

for ith_depth = 1:length(cut_distance)    
    for ith_polytope = 1:length(vertexSkeleton(ith_depth).polytope)
        Nvertices = length(vertexSkeleton(ith_depth).polytope(ith_polytope).vertices(:,1));
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).vector_direction_of_unit_cut(:,1))==Nvertices);
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).unit_vertex_projection_vectors(:,1))==Nvertices);
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).flag_vertexIsNonConvex(:,1))==Nvertices);
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).intersection_points(:,1))==Nvertices);
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).min_cut(:,1))==1);
        assert(length(vertexSkeleton(ith_depth).polytope(ith_polytope).indices_repeated(:,1))==length(vertexSkeleton(ith_depth).polytope(ith_polytope).boundaryEngagedAtMinCut(:,1)));

    end
end

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));


%% Compare speeds of pre-calculation versus post-calculation versus a fast variant
fig_num = 9903;
figure(fig_num);
close(fig_num);

Niterations = 20;

vertices = [0 0; 10 0; 5 10; 0 5; 0 0];

% Do calculation without pre-calculation
tic;
for ith_test = 1:Niterations
    [cut_distance, vertexSkeleton] = fcn_VSkel_polytopeFindVertexSkeleton(vertices,[]);
end
slow_method = toc;

% Do calculation with pre-calculation, FAST_MODE on
tic;
for ith_test = 1:Niterations
    [cut_distance, vertexSkeleton] = fcn_VSkel_polytopeFindVertexSkeleton(vertices,-1);
end
fast_method = toc;

% Plot results as bar chart
figure(373737);
clf;
X = categorical({'Normal mode','Fast mode'});
X = reordercats(X,{'Normal mode','Fast mode'}); % Forces bars to appear in this exact order, not alphabetized
Y = [slow_method fast_method ]*1000/Niterations;
bar(X,Y)
ylabel('Execution time (Milliseconds)')


% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));
