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


%% Demonstration case 1: BASIC example 
fig_num = 1001;
figure(fig_num);
clf;

clear polytopes;
polytopes(1).vertices = [0 0; 4 2; 2 4; -1 3; -2 1; 0 0];

line_spec = '-';
line_width = 3;
fcn_VSkel_plotPolytopes(polytopes,fig_num,line_spec,line_width);


% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Demonstration case 2: BASIC example changing LINE COLOR
fig_num = 1002;
figure(fig_num);
clf;

clear polytopes;
polytopes(1).vertices = [0 0; 4 2; 2 4; -1 3; -2 1; 0 0];

% [FIG]=fcn_VSkel_plotPolytopes(POLYTOPES,FIG_NUM,LINE_SPEC,LINE_WIDTH,COLOR)
% allows the user to specify the input:
% COLOR: a 1-by-3 vector to specify the RGB plot colors [Red Green Blue],
% where 0 <= Red,Green,Blue <= 1

line_spec = '-';
line_width = 3;
line_color1 = [0.5 0 0];
fcn_VSkel_plotPolytopes(polytopes,fig_num,line_spec,line_width,line_color1);


% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Demonstration case 3: BASIC example changing line color, changing LINE WIDTH
fig_num = 1003;
figure(fig_num);
clf;

clear polytopes;
polytopes(1).vertices = [0 0; 4 2; 2 4; -1 3; -2 1; 0 0];

line_spec = '-';
line_width = 7;

fcn_VSkel_plotPolytopes(polytopes,fig_num,line_spec,line_width);


% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));


%% Demonstration case 4: BASIC example changing line color, OPENING A NEW FIGURE BY ENTERING A BLANK FIG_NUM
fig_num = 1004;
figure(fig_num);
clf;

fig_num = []; % Will default to the next figure, usually 1

clear polytopes;
polytopes(1).vertices = [0 0; 4 2; 2 4; -1 3; -2 1; 0 0];

line_spec = '-';
line_width = 3;

% Save all the figures that are currently open
oldFigHandles = get(groot, 'Children');

fig1 = fcn_VSkel_plotPolytopes(polytopes, fig_num, line_spec, line_width);

% Make sure plot opened up a NEW figure
newFigure = gcf;
assert(~any(oldFigHandles==newFigure));


%% Demonstration case 5: BASIC example changing line color, changing AXIS_LIMITS
fig_num = 1005;
figure(fig_num);
clf;

clear polytopes;
polytopes(1).vertices = [0 0; 4 2; 2 4; -1 3; -2 1; 0 0];

line_spec = '-';
line_width = 3;
axis_limits = [0 1 0 1];

fig1 = fcn_VSkel_plotPolytopes(polytopes, fig_num, line_spec, line_width, axis_limits);

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Demonstration case 6: BASIC example changing AXIS_STYLE
fig_num = 1006;
figure(fig_num);
clf;

clear polytopes;
polytopes(1).vertices = [0 0; 4 2; 2 4; -1 3; -2 1; 0 0];

line_spec = '-';
line_width = 3;
line_color1 = [0 1 0]; % Green
axis_limits = [-1 10 -1 10];
axis_style = 'square';

fig1 = fcn_VSkel_plotPolytopes(polytopes, fig_num, line_spec, line_width, axis_limits, axis_style);

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Demonstration case 7: BASIC example changing FILL_INFO
fig_num = 1007;
figure(fig_num);
clf;

clear polytopes;
polytopes(1).vertices = [0 0; 4 2; 2 4; -1 3; -2 1; 0 0];
polytopes(1).cost = 0.4;

line_spec = '-';
line_width = 3;
line_color1 = [0 1 0]; % Green
axis_limits = [-10 10 -10 10];
axis_style = 'square';
fill_info = [1 0 0 1 0.5];

fig1 = fcn_VSkel_plotPolytopes(polytopes, fig_num, line_spec, line_width, axis_limits, axis_style, fill_info);

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

fig_num1 = 9001;
fig_num2 = 9002;
fig_num3 = 9003;
fig_num4 = 9004;
figure(fig_num1); clf;
figure(fig_num2); clf;
figure(fig_num3); clf;
figure(fig_num4); clf;

load('testData_fcn_VSkel_plotPolytopes.mat','polytopes','polytopes2');

line_style = '-';
line_width = 3;
line_color1 = [0.5 0 0];
line_color2 = [0 0 0.5];
line_color3 = [0 0.5 0];
line_color4 = [0 0 0];
axis_box = [0 1 0 1];
axis_style = 'square';
fill_info = [1 0 0 0 0.5];

fig1 = fcn_MapGen_plotPolytopes(polytopes, fig_num1, line_style, line_width, line_color1, axis_box);
fig2 = fcn_MapGen_plotPolytopes(polytopes, fig_num2, line_style, line_width, line_color2, axis_box);
fig3 = fcn_MapGen_plotPolytopes(polytopes, fig_num3, line_style, line_width, line_color3, axis_box, axis_style);
fig4 = fcn_MapGen_plotPolytopes(polytopes, fig_num4, line_style, line_width, line_color4, axis_box, axis_style, fill_info);

% Make sure plots opened up
allFigHandles = get(groot, 'Children');
assert(any(allFigHandles==fig_num1));
assert(any(allFigHandles==fig_num2));
assert(any(allFigHandles==fig_num3));
assert(any(allFigHandles==fig_num4));

% To show that overplotting works, we redo plotting but with a different
% set of polytopes, on the same figures as before
fcn_MapGen_plotPolytopes(polytopes2,fig1,'r-',2);
fcn_MapGen_plotPolytopes(polytopes2,fig2,'b--',2, axis_box);
fcn_MapGen_plotPolytopes(polytopes2,fig3,'g-',3, axis_box,'square');
fcn_MapGen_plotPolytopes(polytopes2,fig4,'k-',3 ,axis_box,'square',[1 0 0 0 0.5]);
