function swain_corban_huang_mod

end



%% Corban Swain Utilities
function fighand = makefigure(fs)
if isfield(fs, 'position')
    fighand = setupfig(fs.n, fs.title, fs.position);
else
    fighand = setupfig(fs.n, fs.title);
end
for i = 1:length(fs.plots)
    subplot(fs.sub(1), fs.sub(2), i);
    hold on;
    makeplot(fs.plots{i});
end
end

function makeplot(ps)
grid on;
if iscell(ps.x) && iscell(ps.y)
    n = length(ps.x);
    cycle = 1;
    spec = {'-'};
    fieldstr = 'spec_cycle';
    if isfield(ps, fieldstr)
        spec = ps.(fieldstr);
    end
    if  n ~= length(ps.y)
        error('x and y cell arrays have mismatched dimensions.');
    end
    for i = 1:n
        fieldstr = 'color_cycle';
        if isfield(ps, fieldstr)
            if (mod(i, ps.(fieldstr)) == 1) && (i ~= 1)
                ax = gca;
                ax.ColorOrderIndex = 1;
                cycle = cycle + 1;
            end
        end
            plot(ps.x{i}, ps.y{i},spec{cycle});
    end
else
    if iscell(ps.x) || iscell(ps.y)
        error('x and y values must both be matrices or cell arrays.');
    end
    plot(ps.x, ps.y);
end
xlabel(ps.xlabel);
ylabel(ps.ylabel);
xlim(ps.xlim);
ylim(ps.ylim);
leg = legend(ps.legend);
if isfield(ps, 'legend_title')
    title(leg, ps.legend_title);
end
leg.LineWidth = 0.5;
if isfield(ps, 'legend_loc')
    leg.Location = ps.legend_loc;
end
end

function new_fig = setupfig(n,title,location)
% SETUPFIGURE Sets up a new figure.
%
% n: figure number
% title: figure title
% location: figure position, [left bottom width height]
%
new_fig = figure(n); clf; hold on;
box off; grid on;
new_fig.Name = title;
if nargin > 2
    new_fig.Position = location;
end
end

function savefig(fig,fig_name)
% SAVEFIGURE Saves the passed figure as a 300 dpi png.

if ~isdir([pwd filesep 'Figures'])
    mkdir 'Figures'
end
f = gobjects(1,1);
name = '';
switch nargin
    case 0
        f = gcf;
    case 1
        f = fig;
    case 2
        f = fig;
        name = fig_name;
end
if isempty(name)
    if isempty(f.Name)
        name = 'Untitled';
    else
        name = fig.Name;
    end
else
    if ~isempty(f.Name)
        name = [name, '-', f.Name];
    end
end
filename = ['Figures' filesep name];
print(f,filename,'-dpng','-r300');
end

function corban_figure_defaults
% CORBANFIGUREDEFAULTS Sets default values to make pretty figures.
fontSize = 13;
font = 'Helvetica';
set(groot, ...
    'defaultLineMarkerSize', 40,...
    'defaultLineLineWidth', 3, ...
    'defaultAxesFontSize', fontSize, ...
    'defaultAxesTitleFontWeight', 'normal', ...
    'defaultAxesFontName', font, ...
    'defaultAxesLabelFontSizeMultiplier', 1.1, ...
    'defaultAxesLineWidth', 2, ...
    'defaultFigureColor', [1 1 1], ...
    'defaultTextInterpreter', 'tex', ...
    'defaultTextFontSize',fontSize, ...
    'defaultTextFontName', font ...
    );
end

function cleanup
clc;
clear;
end