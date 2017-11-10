function swain_corban_huang

function main
    cleanup; % close all;
    corban_figure_defaults;
    figstosim = 4:5;
    for i = figstosim
       fh = makefigure(figures{i}());
       figure(fh)
       savefig(fh);
    end
end

%% Global Variables
global Cf18_source Cfd_source Comfd_source;
omfd_fit_params = [];

%% Importing Data
csv = importdata('plasma_totalF18_FIG3.csv');
Cf18_fig3.t = csv.data(:, 1);
Cf18_fig3.y = csv.data(:, 2);
csv = importdata('plasma_3OMFDFraction_FIG4.csv');
omfdfrac_fig4.t = csv.data(:, 1);
omfdfrac_fig4.y = csv.data(:, 2);
csv = importdata('plasma_totalF18_FIG5.csv');
Cf18_fig5.t = csv.data(:, 1);
Cf18_fig5.y = csv.data(:, 2);
csv = importdata('tissue_radioactivity_striatum_FIG6.csv');
Cs_fig6.t = csv.data(:, 1);
Cs_fig6.y = csv.data(:, 2);
csv = importdata('tissue_radioactivity_cerebellum_FIG6.csv');
Cc_fig6.t = csv.data(:, 1);
Cc_fig6.y = csv.data(:, 2);

%% Figures

figures{1} = @testfig1;
    function testfig1
        Cf18_source = Cf18_fig3;
        setupfig(2, ' Test');
        t2 = linspace(Cf18_source.t(1),Cf18_source.t(end),10000);
        plot(log(Cf18_source.t),Cf18_source.y,'k.',log(t2),Cf18(t2),'-r');
    end

figures{4} = @fig4;
    function fs =  fig4
        fignum = 4;
        params0 = [0.0110, 0.025, 0.0105];
        Cf18_source = Cf18_fig3;
        function omfdfrac = modelfun(p, t)
            Comfd = modelfun_omfd(p, t, 1);
            omfdfrac = Comfd ./ Cf18(t);
        end
        t = omfdfrac_fig4.t; y = omfdfrac_fig4.y;
        loadin = true;
        varname = 'fit_params'; 
        filename = 'omfd_fit_params_cache.mat';
        if loadin            
            load(filename, varname);            
        else
            fit_params = nlinfit(t, y, @modelfun, params0, getfitopts);
            fit_params = abs(fit_params);
            save(filename, varname);
        end
        clear varname filename;
        omfd_fit_params = fit_params;
        setupfig(fignum, 'Figure 4');
        tfit = omfdfrac_fig4.t + eps;
        yfit = modelfun(fit_params, tfit);
        ps.x = {tfit, t};  ps.y = {yfit, y};
        ps.xlim = [0 120]; ps.ylim = [0 1];
        ps.xlabel = 'Time (min)'; ps.ylabel = '3-OMFD Fraction';
        ps.legend = {'Fitting Result', 'Measured Data'};
        ps.legend_loc = 'east';
        ps.linespec = {'-','ko'};
        ps.markersize{2} = 10; 
        ps.cycle_len = 1;
        fs.n = fignum;
        fs.position = [5 613 519 349];
        fs.title = sprintf('%d - Fitting 3-OMFD Frction Timecourse', ...
            fignum);
        fs.plots = {ps};
    end

figures{5} = @fig5;
    function fs = fig5
        fignum = 5;
        fit_params = omfd_fit_params;
        Cf18_source = Cf18_fig5;
        t = Cf18_fig5.t(2:end);
        Comfd = modelfun_omfd(fit_params, t, 1);
        Cfdopa = Cf18(t) - Comfd;
        ps.x = {Cf18_fig5.t, t, t};
        ps.y = {Cf18_fig5.y, Cfdopa, Comfd};
        ps.xlim = [0 120]; ps.ylim = [0.001 10];
        ps.grid = 'off';
        ps.xlabel = 'Time (min)'; 
        ps.ylabel = 'Radioactivity (\muCi/mL)';
        ps.yscale = 'log';
        ps.legend = {'Total F-18', 'FDOPA', '3-OMFD'};
        ps.cycle_len = 1;
        ps.linespec = {'ks-','-','-'};
        ps.markersize = {10, [], []};
        fs.n = fignum;
        fs.position = [526 612 525 351];
        fs.title = sprintf('%d - Plasma Time-Activity Curves', ...
            fignum);
        fs.plots = {ps};
    end

main;
end

%% Settings
function fitopts = getfitopts
fitopts = statset('RobustWgtFun', 'fair', ...
    'Display', 'final');
end

%% Data Interpolators
function yq = interp_templ(x, y, xq)
yq = interp1(x, y, xq, 'pchip');
end

function cq = Cf18(tq)
global Cf18_source;
t = Cf18_source.t;
c = Cf18_source.y;
cq = interp_templ(t, c, tq);
end

function cq = Cfd(tq)
global Cfd_source;
t = Cfd_source.t;
c = Cfd_source.y;
cq = interp_templ(t, c, tq);
end

function cq = Comfd(tq)
global Comfd_source;
t = Comfd_source.t;
c = Comfd_source.y;
cq = interp_templ(t, c, tq);
end

%% ODE Functions
function dY = odefun_striatum(t, Y, p)

% 1: C1
% 2: C2
% 3: C3

p = num2cell(p);
[K1, k2, k3, k4] = p{:};
K5 = K1 * 1.7;
k6 = K5;
A = [-(k2 + k3),   0,   0;
             k3, -k4,   0;
              0,   0, -k6];
B = [  K1 * Cfd(t);
                 0;
     K5 * Comfd(t)];
dY = A * Y + B;
end

function dY = odefun_cerebellum(t, Y, p)
    
% 1: Cp1
% 2: Cp2

p = num2cell(p);
[Kp1, kp2] = p{:};
Kp5 = Kp1 * 1.7;
kp6 = Kp5;
A = [-kp2,    0;
        0, -kp6];
B = [  Kp1 * Cfd(t);
     Kp5 * Comfd(t)];
dY = A * Y + B;
end

function dY = odefun_omfd(t, Y, p)

% 1: Comfd
% 2: Cx

p = num2cell(p);
[kb12, kb2, kb3] = p{:};
A = [-(kb12 + kb2),  kb3;
               kb2, -kb3];
B = [kb12 * Cf18(t);
                  0];
dY = A * Y + B;
end

%% Model Functions
function opts = getodeopts(nvars)   
tol = 1e-6;
opts = odeset('RelTol', tol,...
    'AbsTol',tol, ...
    'NonNegative', 1:nvars);    
end

function Y = modelfun_template(p, t, y_ind_out, nvars, odefun)
tspan = t;
y0 = zeros(1, nvars);
odeopts = getodeopts(nvars);
[~, y_out] = ode45(odefun, tspan, y0, odeopts, abs(p));
if isempty(y_ind_out); y_ind_out = 1:nvars; end
Y = y_out(:, y_ind_out);
end

function Y = modelfun_striatum(p, t, y_ind_out)
if nargin == 2; y_ind_out = []; end
Y = modelfun_template(p, t, y_ind_out, 3, @odefun_striatum);
end

function Y = modelfun_cerebellum(p, t, y_ind_out)
if nargin == 2; y_ind_out = []; end
Y = modelfun_template(p, t, y_ind_out, 2, @odefun_cerebellum);
end

function Y = modelfun_omfd(p, t, y_ind_out)
if nargin == 2; y_ind_out = []; end
Y = modelfun_template(p, t, y_ind_out, 2, @odefun_omfd);
end

%% Corban Swain Utilities
function ret = NIL
ret = "**NIL**";
end

function bool = isnil(x)    
bool = (x == NIL);
end

function [val, S] = qgetfield(S, fieldname, default)
if ~isfield(S, fieldname)
    S.(fieldname) = default;
end
val = S.(fieldname);
end

function fighand = makefigure(fs)
    function val = qgf(fieldname, default)
        [val, fs] = qgetfield(fs, fieldname, default);
    end
fighand = setupfig(fs.n, fs.title, qgf('position', 'auto'));
qgf('sub', [1, 1]);
if prod(fs.sub) < length(fs.plots)
    error('The subplot dimensions are too small for the number of plots.');
end
for i = 1:length(fs.plots)
    subplot(fs.sub(1), fs.sub(2), i);
    hold on;
    makeplot(fs.plots{i});
end
end

function makeplot(ps)
    function val = qgf(fieldname, default)
        [val, ps] = qgetfield(ps, fieldname, default);
    end
isf = @(fieldname) isfield(ps, fieldname);
if ~all(isfield(ps, {'x', 'y'}))
   error('X and Y values must be passed in the plotspec struct.'); 
end
if iscell(ps.x) && iscell(ps.y)
    n = length(ps.x);
    cycle = 1;
    qgf('linespec', {'-'});
    if  n ~= length(ps.y)
        error('x and y cell arrays have mismatched dimensions.');
    end
    for i = 1:n
        fieldname = 'color_order_len';
        if isf(fieldname)
            if (mod(i, ps.(fieldname)) == 1) && (i ~= 1)
                ax = gca;
                ax.ColorOrderIndex = 1;
            end
        end
        fieldname = 'cycle_len';
        if isf(fieldname)
            if (mod((i - 1), ps.(fieldname)) == 0) && (i ~= 1)
                cycle = cycle + 1;
            end
        end             
        % FIXME - Implement linespec like markersize
        p = plot(ps.x{i}, ps.y{i}, ps.linespec{cycle});
        fieldname = 'markersize';
        if isf(fieldname)
            sizespec = ps.(fieldname){i};
            if ~isempty(sizespec)
                p.MarkerSize = sizespec;
            end
        end
            
    end
else
    if iscell(ps.x) || iscell(ps.y)
        error('x and y values must both be matrices or cell arrays.');
    end
    % FIXME - Check for Linespec as cell array ... actually use single
    % value above under all conditions if a cell array is not passed. ...
    % maybe convert x and y to a 1 X 1 cell array to keep everything in the
    % same place...
    plot(ps.x, ps.y, qgf('linespec', '-'));
end
ax = p.Parent;
grid(ax, qgf('grid', 'off'));
ax.XScale = qgf('xscale', 'linear');
ax.YScale = qgf('yscale', 'linear');
xlabel(qgf('xlabel','x'));
ylabel(qgf('ylabel','y'));
xlim(qgf('xlim','auto'));
ylim(qgf('ylim','auto'));
fieldname = 'legend';
if isf(fieldname')
    leg = legend(ps.(fieldname));
    fieldname = 'legend_title';
    if isf(fieldname')
        title(leg, ps.(fieldname));
    end
    leg.LineWidth = qgf('legend_linewidth', 0.5);
    leg.Location = qgf('legend_loc','best');
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
    if location ~= 'auto'
        new_fig.Position = location;
    end
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