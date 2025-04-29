% script_test_fcn_VSkel_polytopeFindVertexAngles
% Tests function: fcn_VSkel_polytopeFindVertexAngles

% REVISION HISTORY:
% 2021_08_01
% -- first written by S. Brennan 
% 2025_04_29 by Sean Brennan
% -- copied test script out of MapGen and into VSkel
%%%%%%%%%%%%%%§

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


%% Basic example of vertex calculation - a square
fig_num = 9001;
figure(fig_num);
clf;

vertices = [0 0; 1 0; 1 1; 0 1; 0 0];
angles =...
    fcn_VSkel_polytopeFindVertexAngles(...
    vertices,fig_num);

assert(1000*eps>abs(360-sum(angles)*180/pi));
% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Basic example of vertex calculation - a square, NO FIGURE
fig_num = 9002;
figure(fig_num);
close(fig_num);

vertices = [0 0; 1 0; 1 1; 0 1; 0 0];
angles =...
    fcn_VSkel_polytopeFindVertexAngles(...
    vertices,[]);

assert(1000*eps>abs(360-sum(angles)*180/pi));
% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));

%% Basic example of vertex calculation - a square, NO FIGURE, FAST MODE
fig_num = 9003;
figure(fig_num);
close(fig_num);

vertices = [0 0; 1 0; 1 1; 0 1; 0 0];
angles =...
    fcn_VSkel_polytopeFindVertexAngles(...
    vertices,-1);

assert(1000*eps>abs(360-sum(angles)*180/pi));
% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));

%% Basic example of vertex calculation - a triangle
fig_num = 9004;
figure(fig_num);
clf;

vertices = [0 0; 1 1; 0 1; 0 0];
angles =...
    fcn_VSkel_polytopeFindVertexAngles(...
    vertices,fig_num);

assert(1000*eps>abs(360-sum(angles)*180/pi));
% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Random polytope calculation
fig_num = 9005;
figure(fig_num);
clf;

vertices = [0 0; 2 2; 0 1; 0 0];

% Basic example of vertex calculation
angles =...
    fcn_VSkel_polytopeFindVertexAngles(...
    vertices,fig_num);
assert(1000*eps>abs(360-sum(angles)*180/pi));
% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));


  