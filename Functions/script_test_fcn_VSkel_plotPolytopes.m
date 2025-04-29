% script_test_fcn_VSkel_plotPolytopes
% Tests function: fcn_VSkel_plotPolytopes

% REVISION HISTORY:
% 2021_06_07
% -- first written by S. Brennan.
% 2023_02_22
% -- better examples
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

URHERE

%% Demonstration case 1: BASIC example 
fig_num = 9001;
figure(fig_num);
clf;

% [FIG_HANDLE]=fcn_VSkel_plotPolytopes(POLYTOPES,FIG_NUM,LINE_SPEC,LINE_WIDTH)
% returns:
% a figure with the polytopes plotted as specified
%
% with inputs:
% POLYTOPES: a 1-by-n seven field structure, where n <= number of polytopes
%   with fields:
%   vertices: a m+1-by-2 matrix of xy points with row1 = rowm+1, where m is
%     the number of the individual polytope vertices
%   xv: a 1-by-m vector of vertice x-coordinates
%   yv: a 1-by-m vector of vertice y-coordinates
%   distances: a 1-by-m vector of perimeter distances from one point to the
%     next point, distances(i) = distance from vertices(i) to vertices(i+1)
%   mean: average xy coordinate of the polytope
%   area: area of the polytope
%   max_radius: distance from the mean to the farthest vertex
% FIG_NUM: figure number to plot the values on
% LINE_SPEC: a string, line specifications such as color, 'r', and line or
% point type, '--'
% LINE_WIDTH: width of the line to be plotted
% By default the axes will be determined by the plot function

clear polytopes;

polytopes(1).vertices = [0 0; 4 2; 2 4; -1 3; -2 1; 0 0];
fig_num = 1;
line_spec = '-';
line_width = 3;
fcn_VSkel_plotPolytopes(polytopes,fig_num,line_spec,line_width);


% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Demonstration case 2: BASIC example changing line color
fig_num = 9001;
figure(fig_num);
clf;


% [FIG]=fcn_VSkel_plotPolytopes(POLYTOPES,FIG_NUM,LINE_SPEC,LINE_WIDTH,COLOR)
% allows the user to specify the input:
% COLOR: a 1-by-3 vector to specify the RGB plot colors [Red Green Blue],
% where 0 <= Red,Green,Blue <= 1


fig_num = 2;
line_color = [0.5 0 0];
fcn_VSkel_plotPolytopes(polytopes,fig_num,line_spec,line_width,line_color);


% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Advanced examples
line_width = 3;
axis_box = [0 1 0 1];
halton_range = [1 100]; % Sets the indices to use from halton set
[polytopes] = fcn_VSkel_haltonVoronoiTiling(halton_range);


% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));


% [FIG]=fcn_VSkel_plotPolytopes(POLYTOPES,FIG_NUM,LINE_SPEC,LINE_WIDTH,AXIS_LIMITS)
% allows the user to specify the input:
% AXIS_LIMTS: a 1-by-4 vector to specify the start and end of the x and y
% axes, [xstart xend ystart yend]
fig_num = []; % Will default to the next figure, usually 1
fig1 = fcn_VSkel_plotPolytopes(polytopes,fig_num,'-',line_width,[0.5 0 0]);


% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%
% [FIG]=fcn_VSkel_plotPolytopes(POLYTOPES,FIG_NUM,LINE_SPEC,LINE_WIDTH,AXIS_STYLE)
% allows the user to specify the input:
% AXIS_STYLE: controls the style of the axis, such as square or equal
fig_num = 998;
fig2 = fcn_VSkel_plotPolytopes(polytopes,fig_num,'-',line_width,[0 0 0.5],axis_box);


% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

% [FIG]=fcn_VSkel_plotPolytopes(POLYTOPES,FIG_NUM,LINE_SPEC,LINE_WIDTH,FILL_INFO)
% allows the user to specify the input:
% FILL_INFO: a 1-by-5 vector to specify wether or not there is fill, the
% color of fill, and the opacity of the fill [Y/N, R, G, B, alpha]
fig_num = 999;
fig3 = fcn_VSkel_plotPolytopes(polytopes,fig_num,'-',line_width,[0 0.5 0],axis_box,'square');

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

% [FIG]=fcn_VSkel_plotPolytopes(POLYTOPES,FIG_NUM,LINE_SPEC,LINE_WIDTH,COLOR,AXIS_LIMITS,AXIS_STYLE,FILL_INFO)
% allows the user to specify any combination of all four inputs in any
% order after LINE_WIDTH
fig_num = 1000;
fig4 = fcn_VSkel_plotPolytopes(polytopes,fig_num,'-',line_width,[0 0 0],axis_box,'square',[1 0 0 0 0.5]);


% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));



% To show that overplotting works, we redo plotting but with a different
% set of polytopes, on the same figures as before
halton_range = [101 200]; % Sets the indices to use from halton set
[polytopes2] = fcn_VSkel_haltonVoronoiTiling(halton_range);
fcn_VSkel_plotPolytopes(polytopes2,fig1,'r-',2);
fcn_VSkel_plotPolytopes(polytopes2,fig2,'b--',2,axis_box);
fcn_VSkel_plotPolytopes(polytopes2,fig3,'g-',3,axis_box,'square');
fcn_VSkel_plotPolytopes(polytopes2,fig4,'k-',3,axis_box,'square',[1 0 0 0 0.5]);
