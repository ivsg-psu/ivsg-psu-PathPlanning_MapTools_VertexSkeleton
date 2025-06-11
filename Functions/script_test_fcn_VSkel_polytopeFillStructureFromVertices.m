% script_test_fcn_VSkel_polytopeFillStructureFromVertices
% Tests function: fcn_VSkel_polytopeFillStructureFromVertices

% REVISION HISTORY:
% 2025_05_24 by Sean Brennan
% -- first written by S. Brennan

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

%% Demonstration case 1: multiple 2D polytope
fig_num = 0001;
figure(fig_num);
clf;

% Call the function
clear vertices
vertices{1} = [0 0; 1 0; 1 1];
vertices{2} = [2 2; 3 4; 1 5];
polytopeStructure = fcn_VSkel_polytopeFillStructureFromVertices(vertices, (fig_num));

% Check variable types
assert(isstruct(polytopeStructure));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Demonstration case 3: merged edges in 2D
fig_num = 0003;
figure(fig_num);
clf;

clear vertices
vertices{1} = [0 0; 1 0; 1 1];
vertices{2} = [0 0; 1 1; 0 1];
polytopeStructure = fcn_VSkel_polytopeFillStructureFromVertices(vertices, (fig_num));

% Check variable types
assert(isstruct(polytopeStructure));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));



%% Basic testing examples in 2D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  ____            _        _______        _   _                ______                           _                       ___  _____
% |  _ \          (_)      |__   __|      | | (_)              |  ____|                         | |                     |__ \|  __ \
% | |_) | __ _ ___ _  ___     | | ___  ___| |_ _ _ __   __ _   | |__  __  ____ _ _ __ ___  _ __ | | ___  ___    ______     ) | |  | |
% |  _ < / _` / __| |/ __|    | |/ _ \/ __| __| | '_ \ / _` |  |  __| \ \/ / _` | '_ ` _ \| '_ \| |/ _ \/ __|  |______|   / /| |  | |
% | |_) | (_| \__ \ | (__     | |  __/\__ \ |_| | | | | (_| |  | |____ >  < (_| | | | | | | |_) | |  __/\__ \            / /_| |__| |
% |____/ \__,_|___/_|\___|    |_|\___||___/\__|_|_| |_|\__, |  |______/_/\_\__,_|_| |_| |_| .__/|_|\___||___/           |____|_____/
%                                                       __/ |                             | |
%                                                      |___/                              |_|
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Basic%20Testing%20%20Examples%20%20-%202D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Basic case 1: single 2D polytope with 3 vertices
fig_num = 1001;
figure(fig_num);
clf;

% Call the function
vertices = [0 0; 1 0; 1 1];
polytopeStructure = fcn_VSkel_polytopeFillStructureFromVertices(vertices, (fig_num));

% Check variable types
assert(isstruct(polytopeStructure));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));



%% Basic case 2: multiple 2D polytope
fig_num = 1002;
figure(fig_num);
clf;

% Call the function
clear vertices
vertices{1} = [0 0; 1 0; 1 1];
vertices{2} = [2 2; 3 4; 1 5];
polytopeStructure = fcn_VSkel_polytopeFillStructureFromVertices(vertices, (fig_num));

% Check variable types
assert(isstruct(polytopeStructure));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));



%% Basic testing examples in 3D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  ____            _        _______        _   _                ______                           _                       ____  _____
% |  _ \          (_)      |__   __|      | | (_)              |  ____|                         | |                     |___ \|  __ \
% | |_) | __ _ ___ _  ___     | | ___  ___| |_ _ _ __   __ _   | |__  __  ____ _ _ __ ___  _ __ | | ___  ___    ______    __) | |  | |
% |  _ < / _` / __| |/ __|    | |/ _ \/ __| __| | '_ \ / _` |  |  __| \ \/ / _` | '_ ` _ \| '_ \| |/ _ \/ __|  |______|  |__ <| |  | |
% | |_) | (_| \__ \ | (__     | |  __/\__ \ |_| | | | | (_| |  | |____ >  < (_| | | | | | | |_) | |  __/\__ \            ___) | |__| |
% |____/ \__,_|___/_|\___|    |_|\___||___/\__|_|_| |_|\__, |  |______/_/\_\__,_|_| |_| |_| .__/|_|\___||___/           |____/|_____/
%                                                       __/ |                             | |
%                                                      |___/                              |_|
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Basic%20Testing%20%20Examples%20%20-%203D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Works in 3D also
% clear poly
% poly.Vertices = [0 0 0; 1 0 1; 1 1 0; 2 2 0; 3 2 0; 4 3 5];
% poly.Faces = [1 2; 2 3; 3 1; 4 5; 5 6; 6 4];
% poly.FaceVertexCData = [0; 1; 0.5];
% poly.FaceColor = 'none'; % flat
% poly.EdgeColor = 'flat';
% vertexNumbering = (1:length(poly.Vertices(:,1)))';
% colorIndices = mod(vertexNumbering-1,length(colorOrdering(:,1)))+1;
% poly.FaceVertexCData = colorOrdering(colorIndices,:);
% poly.LineWidth = 2;
% patch(poly);

fig_num = 3001;
figure(fig_num);
clf;

% Call the function
clear vertices
vertices{1} = [0 0 0; 1 0 1; 1 1 0];
vertices{2} = [2 2 0; 3 2 0; 4 3 5];
polytopeStructure = fcn_VSkel_polytopeFillStructureFromVertices(vertices, (fig_num));

% Check variable types
assert(isstruct(polytopeStructure));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));




%% Demonstration case: 3D polytope (cube)
fig_num = 2002;
figure(fig_num);
close(fig_num);

% The following is motivated by a square cube. It has 6
% external faces, and each vertex connects 3 faces.
clear vertices
vertices{1} = flipud([0 0 0; 0 1 0; 1 1 0; 1 0 0])*5;  % Bottom XY face
vertices{2} = flipud([0 0 0; 0 0 1; 0 1 1; 0 1 0])*5;  % Bottom YZ face
vertices{3} = flipud([0 0 0; 1 0 0; 1 0 1; 0 0 1])*5;  % Bottom XZ face
vertices{4} = flipud([1 1 1; 0 1 1; 0 0 1; 1 0 1])*5;  % Top XY face
vertices{5} = flipud([1 1 1; 1 0 1; 1 0 0; 1 1 0])*5;  % Top YZ face
vertices{6} = flipud([1 1 1; 1 1 0; 0 1 0; 0 1 1])*5;  % Top XZ face

% Call the function
polytopeStructure = fcn_VSkel_polytopeFillStructureFromVertices(vertices, (fig_num));

% Check variable types
assert(isstruct(polytopeStructure));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));




%% Demonstration case: 3D polytope (square bottom pyramid)
%                    V5
%                   /=\\
%                  /===\ \
%                 /=====\' \
%                /=======\'' \
%               /=========\ ' '\
%              /===========\''   \
%             /=============\ ' '  \
%            /===============\  F3'  \
%           /=======F2 =======\' ' ' ' \
%          /===================\' ' '  ' \
%         /=====================\' '   ' ' V3
%        /=======================\  '   ' /
%       /=========================\   ' /
%      /===========================\'  /
%     V1============================V2
%    

fig_num = 2002;
figure(fig_num);
close(fig_num);

clear vertices
vertices{1} = flipud([0 0 0; 0 2 0; 2 2 0; 2 0 0])*5;  % Bottom XY face
vertices{2} = flipud([0 0 0; 2 0 0; 1 1 2])*5;         % front face
vertices{3} = flipud([2 0 0; 2 2 0; 1 1 2])*5;         % right face
vertices{4} = flipud([2 2 0; 0 2 0; 1 1 2])*5;         % back face
vertices{5} = flipud([0 2 0; 0 0 0; 1 1 2])*5;         % left face

% Call the function
polytopeStructure = fcn_VSkel_polytopeFillStructureFromVertices(vertices, (fig_num));

% Check variable types
assert(isstruct(polytopeStructure));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));


%% Demonstration case: 3D polytope (triangular base pyramid)
%                  3
%                 / \ >
%                /   \  >
%               / (3) \   > 
%              /    1  \ 2 < > 4 
%             /      <  \   /   
%            /   <  (4)  \ /   
%           1-<-----------2
%
%
% Correct answer should be:
% 1     2     3     4
     
fig_num = 2003;
figure(fig_num);
close(fig_num);

clear vertices
vertices{1} = flipud([0 0 0; 1 1 0; 2 0 0])*5;    % Bottom XY face
vertices{2} = flipud([0 0 0; 2 0 0; 1 0.5 2])*5;  % front face
vertices{3} = flipud([2 0 0; 1 1 0; 1 0.5 2])*5;  % right face
vertices{4} = flipud([1 1 0; 0 0 0; 1 0.5 2])*5;  % left face

% Call the function
polytopeStructure = fcn_VSkel_polytopeFillStructureFromVertices(vertices, (fig_num));

% Check variable types
assert(isstruct(polytopeStructure));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Demonstration case: 3D polytope (cube created by adjacent faces)
fig_num = 2004;
figure(fig_num);
close(fig_num);

% The following is motivated by a square cube. It has 6
% external faces, and each vertex connects 3 faces.
clear vertices
vertices{1} = flipud([0 0 0; 0 1 0; 1 1 0])*5;  % Bottom XY face (1 of 2)
vertices{2} = flipud([0 0 0; 1 1 0; 1 0 0])*5;  % Bottom XY face (2 of 2)
vertices{3} = flipud([0 0 0; 0 0 1; 0 1 1; 0 1 0])*5;  % Bottom YZ face
vertices{4} = flipud([0 0 0; 1 0 0; 1 0 1; 0 0 1])*5;  % Bottom XZ face
vertices{5} = flipud([1 1 1; 0 1 1; 0 0 1; 1 0 1])*5;  % Top XY face
vertices{6} = flipud([1 1 1; 1 0 1; 1 0 0; 1 1 0])*5;  % Top YZ face
vertices{7} = flipud([1 1 1; 1 1 0; 0 1 0; 0 1 1])*5;  % Top XZ face

% Call the function
polytopeStructure = fcn_VSkel_polytopeFillStructureFromVertices(vertices, (fig_num));

% Check variable types
assert(isstruct(polytopeStructure));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Demonstration case: 3D polytope (square bottom pyramid cut into cube, requiring point insertion)
%                    V5
%                   /=\\
%                  /===\ \
%                 /=====\' \
%                /=======\'' \
%               /=========\ ' '\
%              /===========\''   \
%             /=============\ ' '  \
%            /===============\  F3'  \
%           /=======F2 =======\' ' ' ' \
%          /===================\' ' '  ' \
%         /=====================\' '   ' ' V3
%        /=======================\  '   ' /
%       /=========================\   ' /
%      /===========================\'  /
%     V1============================V2
%    

fig_num = 2005;
figure(fig_num);
close(fig_num);

clear vertices
% vertices{1} = flipud([0 0 0; 0 2 0; 2 2 0; 2 0 0])*5;  % Pyramid bottom XY face
vertices{1} = ([0 0 0; 2 0 0; 1 1 2])*5;         % Pyramid front face
vertices{2} = ([2 0 0; 2 2 0; 1 1 2])*5;         % Pyramid right face
vertices{3} = ([2 2 0; 0 2 0; 1 1 2])*5;         % Pyramid back face
vertices{4} = ([0 2 0; 0 0 0; 1 1 2])*5;         % Pyramid left face

vertices{5} = flipud([0 0 0; 0 0 1; 0 1 1; 0 1 0])*10;  % Cube bottom YZ face
vertices{6} = flipud([0 0 0; 1 0 0; 1 0 1; 0 0 1])*10;  % Cube bottom XZ face
vertices{7} = flipud([1 1 1; 0 1 1; 0 0 1; 1 0 1])*10;  % Cube top XY face
vertices{8} = flipud([1 1 1; 1 0 1; 1 0 0; 1 1 0])*10;  % Cube top YZ face
vertices{9} = flipud([1 1 1; 1 1 0; 0 1 0; 0 1 1])*10;  % Cube top XZ face

% Call the function
polytopeStructure = fcn_VSkel_polytopeFillStructureFromVertices(vertices, (fig_num));

% Check variable types
assert(isstruct(polytopeStructure));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));


% Check that all unit vectors are unit length
unit_normal_vectors_length = sum(unit_normal_vectors.^2,2).^0.5;
assert(all(abs(unit_normal_vectors_length - 1)<1E-10)==true);

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

%% Basic example - NO FIGURE
fig_num = 9901;
figure(fig_num);
close(fig_num);

% Call the function
clear vertices
vertices{1} = [0 0; 1 0; 1 1];
vertices{2} = [2 2; 3 4; 1 5];
polytopeStructure = fcn_VSkel_polytopeFillStructureFromVertices(vertices, ([]));

% Check variable types
assert(isstruct(polytopeStructure));

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));

%% Basic example of vertex calculation - non-normal wall shrinking, NO FIGURE, FAST MODE
fig_num = 9902;
figure(fig_num);
close(fig_num);


% Call the function
clear vertices
vertices{1} = [0 0; 1 0; 1 1];
vertices{2} = [2 2; 3 4; 1 5];
polytopeStructure = fcn_VSkel_polytopeFillStructureFromVertices(vertices, (-1));

% Check variable types
assert(isstruct(polytopeStructure));

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));

%% Compare speeds of pre-calculation versus post-calculation versus a fast variant
fig_num = 9903;
figure(fig_num);
close(fig_num);

Niterations = 100;


% Call the function
clear vertices
vertices{1} = [0 0; 1 0; 1 1];
vertices{2} = [2 2; 3 4; 1 5];

% Do calculation without pre-calculation
tic;
for ith_test = 1:Niterations
    % Call the function
    polytopeStructure = fcn_VSkel_polytopeFillStructureFromVertices(vertices, ([]));
end
slow_method = toc;

% Do calculation with pre-calculation, FAST_MODE on
tic;
for ith_test = 1:Niterations
    % Call the function
    polytopeStructure = fcn_VSkel_polytopeFillStructureFromVertices(vertices, (-1));
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
