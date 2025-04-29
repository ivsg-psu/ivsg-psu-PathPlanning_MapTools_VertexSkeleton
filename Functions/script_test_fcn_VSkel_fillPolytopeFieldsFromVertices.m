% script_test_fcn_VSkel_fillPolytopeFieldsFromVerticies
% Tests: fcn_VSkel_fillPolytopeFieldsFromVerticies

%
% REVISION HISTORY:
%
% 2021_07_02 by Sean Brennan
% -- first write of script
% 2023_01_15 by Sean Brennan
% -- added figure number
% 2025_04_29 by Sean Brennan
% -- copied test script out of MapGen and into VSkel
%%%%%%%%%%%%%%§

close all;

URHERE

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

%% Demonstration case 1: 
fig_num = 9001;
figure(fig_num);
clf;

line_width = 3;
clear polytopes
polytopes(1).vertices = [0 0; 4 2; 2 4; 0 0];
polytopes = fcn_VSkel_fillPolytopeFieldsFromVertices(polytopes);
fcn_VSkel_plotPolytopes(polytopes,fig_num,'r-',line_width);

assert(isequal(polytopes.vertices,[0,0;4,2;2,4;0,0]));
assert(isequal(polytopes.xv,[0,4,2]));
assert(isequal(polytopes.yv,[0,2,4]));
assert(isequal(round(polytopes.distances,4),[4.4721;2.8284;4.4721]));
assert(isequal(polytopes.mean,[2,2]));
assert(isequal(polytopes.area,6));
assert(isequal(round(polytopes.max_radius,4),2.8284));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Demonstration case 2:
fig_num = 9001;
figure(fig_num);
clf;

polytopes(2).vertices = [10 10; 14 21; 12 41; 10 10];
polytopes = fcn_VSkel_fillPolytopeFieldsFromVertices(polytopes);
fcn_VSkel_plotPolytopes(polytopes,fig_num,'r-',line_width);
