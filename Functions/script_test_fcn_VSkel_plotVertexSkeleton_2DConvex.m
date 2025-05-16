% script_test_fcn_VSkel_plotVertexSkeleton_2DConvex
% Tests function: fcn_VSkel_plotVertexSkeleton_2DConvex

% REVISION HISTORY:
% 2025_04_29 by Sean Brennan
% -- first written by S. Brennan by pulling function out of MapGen_polytopeFindVertexSkeleton
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
% fcn_VSkel_plotVertexSkeleton_2DConvex(vertices, projection_vectors, cut_distance, (fig_num))

%% Demonstration case 1: BASIC example 
fig_num = 1001;
figure(fig_num);
clf;

clear verticies projection_vectors cut_distance
verticies{1} = [
     0     0
    10     0
     5    10
     0     5
     0     0];
verticies{2} = [
    4.2794    3.5355
    4.0314    4.0314
    3.5355    3.5355
    4.2794    3.5355];
verticies{3} = [
    3.9809    3.7200
    3.9809    3.7200];
projection_vectors{1} = [
    1.0000    1.0000
    -1.6180    1.0000
    -0.2740   -1.6882
    1.0000   -0.4142
    1.0000    1.0000];
projection_vectors{2} = [
    -1.6180    1.0000
    -0.2740   -1.6882
    2.4142    1.0000
    -1.6180    1.0000];
projection_vectors{3} = [
    0     0
    0     0];

cut_distance = [
    0
    3.5355
    3.7200];


fcn_VSkel_plotVertexSkeleton_2DConvex(verticies, projection_vectors, cut_distance, (fig_num));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));



%% Advanced examples (9XXX figures)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%              _                               _     ______                           _
%     /\      | |                             | |   |  ____|                         | |
%    /  \   __| |_   ____ _ _ __   ___ ___  __| |   | |__  __  ____ _ _ __ ___  _ __ | | ___  ___
%   / /\ \ / _` \ \ / / _` | '_ \ / __/ _ \/ _` |   |  __| \ \/ / _` | '_ ` _ \| '_ \| |/ _ \/ __|
%  / ____ \ (_| |\ V / (_| | | | | (_|  __/ (_| |   | |____ >  < (_| | | | | | | |_) | |  __/\__ \
% /_/    \_\__,_| \_/ \__,_|_| |_|\___\___|\__,_|   |______/_/\_\__,_|_| |_| |_| .__/|_|\___||___/
%                                                                              | |
%                                                                              |_|
% http://patorjk.com/software/taag/#p=display&f=Big&t=Advanced%20%20%20Examples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

