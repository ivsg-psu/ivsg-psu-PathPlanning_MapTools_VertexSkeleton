% script_test_fcn_VSkel_polytopeCentroidAndArea
% Tests: fcn_VSkel_polytopeCentroidAndArea

%
% REVISION HISTORY:
%
% 2021_07_02 by Sean Brennan
% -- first write of script
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

%% Demonstration case 1: Basic example
fig_num = 9001;
figure(fig_num);
clf;


x = [3; 4; 2; -1; -2; -3; -4; -2; 1; 2; 3];
y = [1; 2; 2; 3; 2; -1; -2; -3; -3; -2; 1];
[Centroid,Area] = fcn_VSkel_polytopeCentroidAndArea([x,y],fig_num);

assert(isequal(round(Centroid,4),[-0.1462,-0.2222]));
assert(isequal(round(Area,4),28.5));

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

if 1==1
    plot(x,y,'g-','linewidth',2)
    hold on
    plot(Centroid(:,1),Centroid(:,2),'kx','linewidth',1)
end

%% Demonstration case 2: Basic example, NO FIGURE
fig_num = 9002;
figure(fig_num);
close(fig_num);


x = [3; 4; 2; -1; -2; -3; -4; -2; 1; 2; 3];
y = [1; 2; 2; 3; 2; -1; -2; -3; -3; -2; 1];
[Centroid,Area] = fcn_VSkel_polytopeCentroidAndArea([x,y],[]);

assert(isequal(round(Centroid,4),[-0.1462,-0.2222]));
assert(isequal(round(Area,4),28.5));

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));

if 1==1
    plot(x,y,'g-','linewidth',2)
    hold on
    plot(Centroid(:,1),Centroid(:,2),'kx','linewidth',1)
end

%% Demonstration case 2: Basic example, NO FIGURE, FAST MODE
fig_num = 9003;
figure(fig_num);
close(fig_num);


x = [3; 4; 2; -1; -2; -3; -4; -2; 1; 2; 3];
y = [1; 2; 2; 3; 2; -1; -2; -3; -3; -2; 1];
[Centroid,Area] = fcn_VSkel_polytopeCentroidAndArea([x,y],-1);

assert(isequal(round(Centroid,4),[-0.1462,-0.2222]));
assert(isequal(round(Area,4),28.5));

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));

if 1==1
    plot(x,y,'g-','linewidth',2)
    hold on
    plot(Centroid(:,1),Centroid(:,2),'kx','linewidth',1)
end