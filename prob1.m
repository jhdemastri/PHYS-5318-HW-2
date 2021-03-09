%Author: John DeMastri
%Date: 3/9/2021
%For: PHYS 5318 HW #2 Question 1

clear; close all; clc;
%% Setup
%Let's instantiate the data
rawdata = readtable("problem1data.txt");
rawxdata = removevars(rawdata, {'y1', 'y2', 'y3', 'y4', 'y5', 'y6', 'y7', 'y8', 'y9', 'y10'});
rawydata = removevars(rawdata, {'x'});

%model, f = Ax^b
fit = @(p, x) p(1) .* x .^ p(2);
%initialize our parameters for the fminsearch later
A    = 1;        % Free parameter #1
b    = 2;        % Free parameter #2
%Domain
x = rawdata.x;

%some constants we know from the data
pts  = 15;       % Data points 
dof  = pts - 2;  % Degrees of freedom of the model (Assume 2 fit params)
Nrep = 10;       % Number of measurements (need at least 2)


%convert y-values to an array, can't do math with a table
yarray = table2array(rawydata);
%find the mean and standard error of the mean for each point
ymean = mean(yarray, 2);
yse = (std(yarray')/sqrt(Nrep))';

%concatenate into a table for plotting later
condata = addvars(rawxdata, ymean, yse);

%% Chisq minimization
%weighted chi squared
chi2 = @(p) sum( ((ymean - fit(p, x)) ./ yse).^2 );

% Initial condition for fminsearch 
x0 = [A b];
 
% Use 'fminsearch' to find the parameters that minimize the chi2 function
pmin = fminsearch(chi2, x0);
 
% Find minimum chi2 using parameters, as well as reduced chi2
minchisq = chi2(pmin);
redchisq = minchisq/dof;
 
% Calculate the "goodness of fit" R^2
Rsq_chi2 = 1 - sum((ymean - fit(pmin, x)).^2)/sum((ymean - mean(ymean)).^2);

%% Chisq + 1
% Take a look at chi-squared in parameter space
% NOTE: Adjust the range 'r' to zoom in/out on the chi^2 + 1 ellipse;
% It should fill at least 1% of the grid. 
r = 0.2;
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
[idxA, idxb] = find(Z <= minchisq + 1);
 
% Make sure range was set correctly by checking how many grid points were
% found in the chi2 + 1 ellipse (percentage of grid points)
if numel(idxA)/numel(X) < 0.01
    warning('Ellipse is too small, please zoom in for better precision.');
end
 
% Compute uncertainty in A and b
dAmin = min(Arange(idxA)); dAmax = max(Arange(idxA));
dbmin = min(brange(idxb)); dbmax = max(brange(idxb));
dA = (dAmax - dAmin)/2;    db = (dbmax - dbmin)/2;

%% OLS and WLS
% linearize the data, in this case we log everything
x_lin = log(x);
y_lin = log(ymean); 
 
% OLS

x_ols = [(x_lin.^0) (x_lin.^1)];
y_ols = y_lin;
 
% param ols = (X'*X)^(-1) * X' * Y;
p_ols = inv(x_ols' * x_ols) * x_ols' * y_ols;
 
% OLS coefficients, corrected for the linearization
A_ols = exp(p_ols(1));
b_ols = p_ols(2);

p_ols = [A_ols b_ols];

% Goodness of fit
ols_fit = fit(p_ols, x);
Rsq_ols = 1 - sum((ymean - ols_fit).^2)/sum((ymean- ones(length(ymean),1)*mean(ymean)).^2);
 
 
% WLS

%correct the uncertainty for the linearization
yse_lin = (yse./ymean).^2;

% x and y get weighted by their uncertainty
x_wls = x_ols ./ yse_lin;
y_wls = y_lin ./ yse_lin;
 
% Covaraince matrix
cov = inv(x_wls'*x_wls);

% Compute coeficients (same as OLS method)
p_wls = cov * x_wls' * y_wls;


% Individual coefficients, corrected for linearization
A_wls = exp(p_wls(1));
b_wls = p_wls(2);

%error on the coefficients
sigma_p = sqrt(diag(cov));
dA_wls = A_wls*sigma_p(1)
db_wls = sigma_p(2)

%parameters, corrected for linearization
p_wls = [A_wls b_wls]
 
% Goodness of fit
wls_fit = fit(p_wls, x);
Rsq_wls = 1 - sum((ymean - wls_fit).^2)/sum((ymean - ones(length(ymean),1)*mean(ymean)).^2);

%% Plot and Output
%make sure we have an appropriately sized figure
figure('Renderer', 'painters', 'Position', [400 300 900 600])

%plot everything on a logscale
errorbar(condata.x, condata.ymean, condata.yse, 'x');
hold on
plot(x, fit(pmin, x), 'k');
plot(x, fit(p_ols, x), '--k');
plot(x, fit(p_wls, x), '-.k');
hold off
set(gca, 'XScale', 'log', 'YScale', 'log')

%labels, legends, and annotations
xlabel('Domain');
ylabel('Range');
title('Toy Model Fit To Dummy Data On a Log Scale');

str_chi2 = sprintf(['\\chi^2 Fit Params:\n'...
               'A = %3.2f \\pm %3.3f \n' ...
               'p = %3.2f \\pm %3.3f \n' ...
               'R^2 = %4.4f \n' ...
               '\\chi^2_{min} = %4.3f \n' ...
               '\\chi^2_{red} = %4.3f \n'], pmin(1), dA, pmin(2),...
                db, Rsq_chi2, minchisq, redchisq);
            
str_ols = sprintf(['OLS Fit Params: \n' ...
               'A = %3.2f, p = %3.2f \n' ...
               'R^2 = %4.4f \n'], A_ols, b_ols, Rsq_ols);
           
str_wls = sprintf(['WLS Fit Params: \n' ...
               'A = %3.2f \\pm %3.3f \n' ...
               'p = %3.2f \\pm %3.3f \n' ...
               'R^2 = %4.4f \n'], A_wls, dA_wls, b_ols, db_wls, Rsq_wls);
legend('Data', sprintf('\\chi^2_{min}'), 'ols', 'wls', 'Location','Northwest');

a = gca; % get the current axis;
% set the width of the axis (the third value in Position) 
% to be 60% of the Figure's width
a.Position(3) = 0.7;            
annotation('textbox', [.84, .75, .1,.1], 'String',  str_chi2); 
annotation('textbox', [.84, .425, .1,.1], 'String',  str_ols);
annotation('textbox', [.84, .225, .1,.1], 'String',  str_wls);
set(gcf, 'Name', 'Problem 1 by John DeMastri', 'NumberTitle', 'Off')