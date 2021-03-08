%Author: John DeMastri
%Date: 3/9/2021
%For: PHYS 5318 HW #2 Question 1

clear; close all; clc;

%% Setup
%Let's instantiate the data
rawdata = readtable('problem2data.txt');
rawtdata = removevars(rawdata, {'x1', 'x2', 'x3', 'Var5'});
rawxdata = removevars(rawdata, {'t', 'Var5'});

%model, b(t) = b_eq - [b_eq -0.355]exp(-k*t)
fit = @(p, t) p(1) - (p(1) - 0.355)*exp(-p(2)* t);

%initialize our parameters for the fminsearch later
b_eq    = 1;        % Free parameter #1
k    = 2;        % Free parameter #2
%Domain
t = rawdata.t;

%some constants we know from the data
pts  = 25;       % Data points 
dof  = pts - 2;  % Degrees of freedom of the model (Assume 2 fit params)
Nrep = 3;       % Number of measurements (need at least 2)


%convert y-values to an array, can't do math with a table
xarray = table2array(rawxdata);
%find the mean and standard error of the mean for each point
xmean = mean(xarray, 2);
xse = (std(xarray')/sqrt(Nrep))';

%concatenate into a table for plotting later
condata = addvars(rawtdata, xmean, xse);

%% Chi Sq Analysis
%weighted chi squared
chi2 = @(p) sum( ((xmean - fit(p, t)) ./ xse).^2 );

% Initial condition for fminsearch 
t0 = [b_eq k];
 
% Use 'fminsearch' to find the parameters that minimize the chi2 function
pmin = fminsearch(chi2, t0);
 
% Find minimum chi2 using parameters, as well as reduced chi2
minchisq = chi2(pmin);
redchisq = minchisq/dof;
 
% Calculate the "goodness of fit" R^2
Rsq_chi2 = 1 - sum((xmean - fit(pmin, t)).^2)/sum((xmean - mean(xmean)).^2);

%% Chisq + 1
% Take a look at chi-squared in parameter space
% NOTE: Adjust the range 'r' to zoom in/out on the chi^2 + 1 ellipse;
% It should fill at least 1% of the grid. 
r = 0.01;
Arange = linspace(pmin(1) - r, pmin(1) + r, 300);
brange = linspace(pmin(2) - r, pmin(2) + r, 300);
[X, Y] = ndgrid(Arange, brange);
 
% Compute Chi-squared at each point in parameter space
Z = zeros(size(X));
for i = 1:size(X, 1)
    for j = 1:size(X, 2)
        Z(i, j) = chi2( [X(i, j), Y(i, j)] );
    end
end
 
% Look for all values in Z that compute as less than Chi^2 + 1
[idxb_eq, idxk] = find(Z <= minchisq + 1);
 
% Make sure range was set correctly by checking how many grid points were
% found in the chi2 + 1 ellipse (percentage of grid points)
if numel(idxb_eq)/numel(X) < 0.01
    warning('Ellipse is too small, please zoom in for better precision.');
end
 
% Compute uncertainty in A and b
db_eq_min = min(Arange(idxb_eq)); db_eq_max = max(Arange(idxb_eq));
dkmin = min(brange(idxk)); dkmax = max(brange(idxk));
db_eq = (db_eq_max - db_eq_min)/2;    dk = (dkmax - dkmin)/2;

%% Plot and Output
%make sure we have an appropriately size
figure('Renderer', 'painters', 'Position', [400 300 900 600])

%plot everything on a logscale
errorbar(condata.t, condata.xmean, condata.xse, 'x');
hold on
plot(t, fit(pmin, t), 'k');
hold off

%labels, legends, and annotations
xlabel('Incubation Time');
ylabel('Extension per Base Pair');
title('Toy Model Fit To Dummy DNA Data');

str_chi2 = sprintf(['\\chi^2 Fit Params:\n' ...
               'b_{eq} = %3.2f \\pm %3.4f \n' ...
               'k = %3.2f \\pm %3.4f \n' ...
               'R^2 = %4.3f \n' ...
               '\\chi^2_{min} = %4.3f \n' ...
               '\\chi^2_{red} = %4.3f \n'], pmin(1), db_eq, pmin(2),...
                dk, Rsq_chi2, minchisq, redchisq);
            
legend('Data', sprintf('\\chi^2_{min}'), 'Location','Northwest');

a = gca; % get the current axis;
% set the width of the axis (the third value in Position) 
% to be 60% of the Figure's width
a.Position(3) = 0.7;            
annotation('textbox', [.835, .75, .1,.1], 'String',  str_chi2); 
set(gcf, 'Name', 'Problem 2 by John DeMastri', 'NumberTitle', 'Off')
