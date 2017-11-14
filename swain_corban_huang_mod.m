%% IMPLEMENTATION OF HUANG, 1991
% MATLAB implementation of the multi-compartment model of FDOPA activity in
% the striatum as observed by PET. This model was originally concieved by
% Sung-cheng Huang and coworkes [1]. Through the following code, figues, 4, 
% 5, and 6B along with the underlying model are reproduced.
% 
% by Corban Swain
% November 7 - 14, 2017
% for 20.420 (MIT Class)
%
% 1. Huang, Sung-Cheng, et al. "Kinetics and modeling of L-6-[18F] 
% fluoro-dopa in human positron emission tomographic studies." Journal of 
% Cerebral Blood Flow & Metabolism 11.6 (1991): 898-913.
%%
function swain_corban_huang
function main
    cleanup; close all;
    fprintf('Beginning Script...\n');
    corban_figure_defaults;
    fig4_calcs; fig5_calcs;
    figstosim = 6;
    for i = figstosim
        fprintf('\nRunning Modification Figure %d %s \n', ... 
            i, repmat('-', 1, 40));
        fh = makefigure(figures{i}());
        figure(fh)
        savefig(fh);
    end
    fprintf('Done!\n');
end

%% Global Variables
global Cf18_source Cfd_source Comfd_source;
omfd_fit_params = [];

%% Importing Data
csv = importdata('wholeblood_totalF18_FIG3.csv');
Cf18_wb_fig3.t = csv.data(:, 1);
Cf18_wb_fig3.y = csv.data(:, 2);
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

num = 60;
t_templ = linspace(eps, 10, num)';
t_templ = [t_templ(1:(end - 1)); linspace(10, 120, num + 1)'];
%% Figures

    function fig4_calcs
        Cf18_source = Cf18_fig3;
        t = omfdfrac_fig4.t; y = omfdfrac_fig4.y;
        function omfdfrac = modelfun(p, t)
            Comfd = modelfun_omfd(p, t, 1);
            omfdfrac = Comfd ./ Cf18(t);
        end
        params0 = [0.0110, 0.025, 0.0105];
        loadin = true;
        varname = 'fit_params'; 
        filename = 'figure4_mod_cache.mat';
        if loadin && isfile(filename)
            load(filename, varname);            
        else
            fprintf(['Beginnging Nonlinear Fit for 3-OMFD Reaction ', ...
                'and Transport Coeffecients...\n']);
            fit_params = nlinfit(t, y, @modelfun, params0, getfitopts);
            fit_params = abs(fit_params);
            save(filename, varname);
        end
        clear varname filename;
        fprintf('Fitted 3-OMFD Reaction Model Params [kb12 kb2 kb3]:\n');
        disp(fit_params); 
        omfd_fit_params = fit_params;
    end

    function fig5_calcs
        fit_params = omfd_fit_params;
        Cf18_source = Cf18_fig5;
        t = Cf18_fig5.t;
        Comfd = modelfun_omfd(fit_params, t, 1);
        Cfdopa = Cf18(t) - Comfd;
        Comfd_source = struct('t', t, 'y', Comfd);
        Cfd_source = struct('t', t, 'y', Cfdopa);
    end

figures{6} = @fig6;
    function fs = fig6
        fignum = 6;
        loadin = true;
        % FIXME - Automate Saving / Checkpoints.
        varnames = {'Cbackground', 'Cwb', 'fit_part_coeff'}; 
        filename = 'figure6_mod_cache_1.mat';
        if loadin && isfile(filename)
            for i = 1:length(varnames)
                load(filename, varnames{i});
            end
        else
            % Fitting Partition Coeffecients
            fprintf(['Beginnging Nonlinear Fit for Plasma/Whole ', ...
                'Blood Partition Coeffecients...\n']);
            part_coef0 = [1.01, 1.08];
            Cwb = @(p, t) abs(p(1)) * 0.6 * Cfd(t) + abs(p(2)) * Comfd(t);
            fit_part_coeff = nlinfit(Cf18_wb_fig3.t, Cf18_wb_fig3.y, ...
                Cwb, part_coef0, getfitopts);
            fit_part_coeff = abs(fit_part_coeff);
            Cbackground = @(t) 0.05 * Cwb(fit_part_coeff, t);
            save(filename, varnames{1});
            for i = 2:length(varnames)
                save(filename, varnames{i}, '-append');
            end
        end
        clear filename varnames;
        fprintf('Fitted Partition Coeffeceints [p1 p2]:\n');
        disp(fit_part_coeff);
        function [stri_activity, C] = modelfun_s(p, t)
            C = modelfun_striatum(p, t);
            stri_activity = sum(C(:, 1:3), 2) + Cbackground(t);
        end
        
        Cwb_mod = @(p, t, Cfa) abs(p(1)) * 0.6 * (Cfd(t) + Cfa) ...
            + abs(p(2)) * Comfd(t);
        Cbackground_mod = @(t, Cfa) 0.05 ... 
            * Cwb_mod(fit_part_coeff, t, Cfa);
        function [stri_activity, C] = modelfun_s_mod(p, t)
            C = modelfun_striatum_mod(p, t);
            stri_activity = sum(C(:, 1:3), 2) ...
                + Cbackground_mod(t,C(:,4));
        end
        
        varnames = {'fit_params', 't', 'y', 'modelfun'};
        filename = 'figure6_mod_cache_2.mat';
        if loadin && isfile(filename)
            for i = 1:length(varnames)
                load(filename, varnames{i});
            end
        else
            % Fitting Cerebellal and Striatal Activities
            fprintf(['Beginnging Nonlinear Fits for ', ...
                'Striatal Activities w/ & w/o Plasma Clearance...\n']);
            n = 1;
            t{n} = Cs_fig6.t; y{n} = Cs_fig6.y;
            modelfun{n} = @modelfun_s;
            params0{n} = [0.0313, 0.0413, 0.0411, 0.0043];  
            n = n + 1;
            t{n} = Cs_fig6.t; y{n} = Cs_fig6.y;
            modelfun{n} = @modelfun_s_mod;
            params0{n} = [0.0313, 0.0413, 0.0411, ([0.0043, 0.0043] ./ 2)];
            fit_params = cell(n, 1);
            for i = 1:n
                fit_params{i} = nlinfit(t{i}, y{i}, modelfun{i}, ...
                    params0{i}, getfitopts);
                fit_params{i} = abs(fit_params{i});
            end
            save(filename, varnames{1});
            for i = 2:length(varnames)
                save(filename, varnames{i}, '-append');
            end
        end
        clear filename varnames;
        fprintf('Fitted INITIAL striatal Model Params [k1 k2 k3 k4]:\n');
        disp(fit_params{1});  
        fprintf(['Fitted MODIFIED striatal Model Params ', ...
            '[k1 k2 k3 k4 k7]:\n']);
        disp(fit_params{2}); 
        n = length(fit_params);
        tout = cell(n, 1); 
        ytot_out = tout; ycomp_out = tout;
        for i = 1:n
            % QUESTION - Plots look weird with t_templ ... why?
            tout{i} = t{i} + eps;
            [ytot_out{i},  ~] ...
                = modelfun{i}(fit_params{i}, tout{i});  
            tout{i} = t_templ;
            [~, ycomp_out{i}] ...
                = modelfun{i}(fit_params{i}, tout{i});
        end
        
        % Preparing plot and figure
        n = 1;
        ps.x{n} = t{1}; ps.y{n} = y{1}; % original striatal
        n = n + 1;
        ps.x{n} = t{1} + eps; ps.y{n} = ytot_out{1}; % fitted striatal
        n = n + 1;
        ps.x{n} = t_templ;
        ps.y{n} = sum(ycomp_out{1}(:,[1, 3]), 2); % 3-omfd + fdopa 
        n = n + 1;
        ps.x{n} = t{2} + eps; ps.y{n} = ytot_out{2}; % fitted striatal, mod
        n = n + 1;
        ps.x{n} = t_templ;
        ps.y{n} = sum(ycomp_out{2}(:,[1, 3]), 2); % 3-omfd + fdopa, mod
        
        ps.xlim = [0 120]; ps.ylim = [0 0.3];
        ps.grid = 'off';
        ps.xlabel = 'Time (min)'; 
        ps.ylabel = 'Radioactivity (\muCi/mL)';
        ps.legend = {'Striatal Activity', 'Striatal Fit', ...
            '3-OMFD & FDOPA', 'Plasma Clearance: Striatal Fit', ...
            'Plasma Clerance: 3-OMFD & FDOPA' };
        ps.legend_loc = 'northeastoutside';
        ps.cycle_len = 1;
        ps.linespec =        {'ks', '-k', 'k:', '-r', ':r'};
        ps.markerfacecolor = { 'k',   [],   [],   [],   []};
        ps.markersize =      {   8,   [],   [],   [],   []};
        ps.linewidth =       {   1,   [],   [],   [],   []};
        fs.n = fignum;
        fs.position = [5 107 893 417];
        fs.title = sprintf(['%d - Striatal ',  ...
            'Activity Timecourses with Modification'], ...
            fignum);
        fs.plots = {ps};
    end

main;
end

%% Settings
function fitopts = getfitopts
fitopts = statset('RobustWgtFun', 'cauchy', ...
    'Display', 'iter', ...
    'DerivStep', eps .^ (1 / 3));
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

function dY = odefun_striatum_mod(t, Y, p)

% 1: C1
% 2: C2
% 3: C3
% 4: Cfa <-- Modification

p = num2cell(p);
%                vv Modification 
[K1, k2, k3, k4, k7] = p{:};
K5 = K1 * 1.7;
k6 = K5;
A = [-(k2 + k3),   0,   0,   0;
             k3, -k4,   0,   0;
              0,   0, -k6,   0;
              0, -k4,   0, -k7]; % Modification
B = [  K1 * Cfd(t);
                 0;
     K5 * Comfd(t);
                 0];
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
tol = 1e-9;
opts = odeset('RelTol', tol,...
    'AbsTol',tol, ...
    'NonNegative', 1:nvars);    
end

function [Y, t] = modelfun_template(p, t, y_ind_out, nvars, odefun)
tspan = t;
y0 = zeros(1, nvars);
odeopts = getodeopts(nvars);
[t, y_out] = ode45(odefun, tspan, y0, odeopts, abs(p));
if isempty(y_ind_out); y_ind_out = 1:nvars; end
Y = y_out(:, y_ind_out);
end

function [Y, t] = modelfun_striatum(p, t, y_ind_out)
if nargin == 2; y_ind_out = []; end
[Y, t] = modelfun_template(p, t, y_ind_out, 3, @odefun_striatum);
end

function [Y, t] = modelfun_striatum_mod(p, t, y_ind_out)
if nargin == 2; y_ind_out = []; end
[Y, t] = modelfun_template(p, t, y_ind_out, 4, @odefun_striatum_mod);
end

function [Y, t] = modelfun_cerebellum(p, t, y_ind_out)
if nargin == 2; y_ind_out = []; end
[Y, t] = modelfun_template(p, t, y_ind_out, 2, @odefun_cerebellum);
end

function [Y, t] = modelfun_omfd(p, t, y_ind_out)
if nargin == 2; y_ind_out = []; end
[Y, t] = modelfun_template(p, t, y_ind_out, 2, @odefun_omfd);
end

%% Corban Swain Utilities

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
        % FIXME - maybe functionalize this more? Be able to take in shorter
        % cell arrays then vars.
        p = plot(ps.x{i}, ps.y{i}, ps.linespec{cycle});
        fieldname = 'markersize';
        if isf(fieldname)
            sizespec = ps.(fieldname){i};
            if ~isempty(sizespec)
                p.MarkerSize = sizespec;
            end
        end
        fieldname = 'markerfacecolor';  
        if isf(fieldname)
            colorspec = ps.(fieldname){i};
            if ~isempty(colorspec)
                p.MarkerFaceColor = colorspec;
            end
        end
        
        fieldname = 'linewidth';  
        if isf(fieldname)
            lw = ps.(fieldname){i};
            if ~isempty(lw)
                p.LineWidth = lw;
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