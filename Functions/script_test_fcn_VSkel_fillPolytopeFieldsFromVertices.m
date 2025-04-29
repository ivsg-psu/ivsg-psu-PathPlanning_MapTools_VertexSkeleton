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

%% Basic Examples (1XXX figures)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  ____            _          ______                           _
% |  _ \          (_)        |  ____|                         | |
% | |_) | __ _ ___ _  ___    | |__  __  ____ _ _ __ ___  _ __ | | ___  ___
% |  _ < / _` / __| |/ __|   |  __| \ \/ / _` | '_ ` _ \| '_ \| |/ _ \/ __|
% | |_) | (_| \__ \ | (__    | |____ >  < (_| | | | | | | |_) | |  __/\__ \
% |____/ \__,_|___/_|\___|   |______/_/\_\__,_|_| |_| |_| .__/|_|\___||___/
%                                                       | |
%                                                       |_|
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Basic%20%20%20Examples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% BASIC EXample
fig_num = 1001;
figure(fig_num);
clf;

clear polytopes
polytopes(1).vertices = [0 0; 4 2; 2 4; 0 0];
is_nonconvex = 0;

polytopes = fcn_VSkel_fillPolytopeFieldsFromVertices(polytopes, (is_nonconvex), (fig_num));

assert(isequal(polytopes.vertices,[0,0;4,2;2,4;0,0]));
assert(isequal(polytopes.xv,[0,4,2]));
assert(isequal(polytopes.yv,[0,2,4]));
assert(isequal(round(polytopes.distances,4),[4.4721;2.8284;4.4721]));
assert(isequal(polytopes.mean,[2,2]));
assert(isequal(polytopes.area,6));
assert(isequal(round(polytopes.max_radius,4),2.8284));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));


%% BASIC EXample, NO FIGURE
fig_num = 1002;
figure(fig_num);
close(fig_num);

clear polytopes
polytopes(1).vertices = [0 0; 4 2; 2 4; 0 0];
is_nonconvex = 0;

polytopes = fcn_VSkel_fillPolytopeFieldsFromVertices(polytopes, (is_nonconvex), ([]));

assert(isequal(polytopes.vertices,[0,0;4,2;2,4;0,0]));
assert(isequal(polytopes.xv,[0,4,2]));
assert(isequal(polytopes.yv,[0,2,4]));
assert(isequal(round(polytopes.distances,4),[4.4721;2.8284;4.4721]));
assert(isequal(polytopes.mean,[2,2]));
assert(isequal(polytopes.area,6));
assert(isequal(round(polytopes.max_radius,4),2.8284));


% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));

%% BASIC EXample, NO FIGURE, FAST MODE
fig_num = 1003;
figure(fig_num);
close(fig_num);

clear polytopes
polytopes(1).vertices = [0 0; 4 2; 2 4; 0 0];
is_nonconvex = 0;

polytopes = fcn_VSkel_fillPolytopeFieldsFromVertices(polytopes, (is_nonconvex), (-1));

assert(isequal(polytopes.vertices,[0,0;4,2;2,4;0,0]));
assert(isequal(polytopes.xv,[0,4,2]));
assert(isequal(polytopes.yv,[0,2,4]));
assert(isequal(round(polytopes.distances,4),[4.4721;2.8284;4.4721]));
assert(isequal(polytopes.mean,[2,2]));
assert(isequal(polytopes.area,6));
assert(isequal(round(polytopes.max_radius,4),2.8284));


% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));

%% Basic case: Two polytopes
fig_num = 1004;
figure(fig_num);
clf;


clear polytopes
polytopes(1).vertices = [0 0; 4 2; 2 4; 0 0];
polytopes(2).vertices = [10 10; 14 21; 12 41; 10 10];
is_nonconvex = 0;

polytopes = fcn_VSkel_fillPolytopeFieldsFromVertices(polytopes, (is_nonconvex), (fig_num));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));