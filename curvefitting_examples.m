%%%%
% This matlab file demonstrates how to fit data to a model using
%   - Chi^2 minimization (with Min chi^2 plus one)
%   - Ordinary least squares (OLS)
%   - Weighted least squares (WLS)
%
% This file is hardcoded for a 2-parameter models, but may be generalized
% for any number of parameters.
%
% Three examples are pre-programmed in for you to explore: a sinusoid
% measurement, a linear fit with signal dependent noise, and the example
% from the counterpart excel sheet on blackboard.
%
% Exploration of this file is encouraged. To guide your investigation,
% think about the following questions.
%
% 1. How do the measurement settings (pts, Nrep, nlvl) affect the precision
%    of the fits?
%
% 2. What is more desirable: more data points per measurement or more total
%    measurements?
%
% 3. How useful is the "Goodness of fit" R^2 statistic at comparing each
%    fitting method?
%
% 4. How do the weighting schemes affect the fits in both methods?
%
% 
% By:  Eric Kercher 
% On:  18FEB2019
% For: PHYS5318 - Principals of Experimental Physics
% At:  Northeastern University
%%%%%%n
 
clear; close all; clc;
 
% Seed for random number generatior. Creates the same noise, and therefore
% the same fit, every run. Comment out for different noise on each run.
rng(5318);
 
% Choose example 'linear' or 'sine' to generate a simulated data set,
% 'excel' to display an example from the excel spreadsheet, 
% or choose 'custom' to write or load in your own data
examp = 'linear';
 
% Choose between 'ordinary' and 'weighted' chi-squared methods
method = 'weighted';
 
 
%% Simulate a measured data set 
A    = 1;        % Free parameter #1
b    = 2;        % Free parameter #2
pts  = 30;       % Data points (shoud be much better than Nyquist!!)
dof  = pts - 2;  % Degrees of freedom of the model (Assume 2 fit params)
Nrep = 5;       % Number of measurements (need at least 2)
nlvl = 0.5;      % Noise level (0 < nlvl < 1)
 
% Domain
x = linspace(0, 2*pi, pts)';
 
% Generate data
switch examp
    case 'sine'
        % Sine wave with Amplitude 'A' and frequency 'b'
        y = repmat(A*sin(b*x), [1 Nrep]) + nlvl*randn(length(x), Nrep);
        
        % Define model using anonymous function
        fit = @(p, x) p(1) .* sin(p(2) .* x);
        
        titlestr = 'Simulated Measurement: Y = Asin(bx)';
        
        % Compute average and SEM
        yavg = mean(y, 2);
        ysem = std(y, 0 , 2) ./ sqrt(Nrep);
 
    case 'linear'
        % Linear response with slope 'A' and intercept 'b'
        % Note the signal dependent noise and constant noise floor.
        y = repmat(A*x + b, [1 Nrep]) + nlvl.*x.*randn(length(x), Nrep) ...
            + .1.*randn(length(x), Nrep);
        
        % Define  model
        fit = @(p, x) p(1) .* x + p(2);
        
        titlestr = 'Simulated Measurement: Y = Ax + b';
        
        % Compute average and SEM
        yavg = mean(y, 2);
        ysem = std(y, 0 , 2) ./ sqrt(Nrep);
 
    case 'excel'
        x    = [1 2 3 4 5]';
        yavg = [0.811938401;
                6.96543519;
                8.350586879;
                9.006213575;
                13.64499783];
        ysem = [0.5 1.2 1.5 1.4 0.8]';
        
        A = 3;
        b = -2;
        pts = length(x);
        dof = pts - 2;
        
        % Define  model
        fit = @(p, x) p(1) .* x + p(2);
        
        titlestr = 'Example from Excel Sheet: Y = Ax + b';
        
    case 'custom'
        % Load your own data here... 
%         load('mydata.mat')
        
%         x = ...
%         y = ...
%         A = ...
%         b = ...
%         pts = ...
%         dof = pts - 2;
%         Nrep = ...
 
%         fit = ...
%         titlestr = ...
 
%         yavg = ...
%         ysem = ...
    
    otherwise
        error('Example not recognized');
end
 
 
%% Chi-squared minimization
 
% Define chi-squared error metric
switch method
    case 'weighted'
        % Data points are weighted by uncertainty
        chi2 = @(p) sum( ((yavg - fit(p, x)) ./ ysem).^2 );
        
    case 'ordinary'
        % Assume equal uncertainties of 1 (This is just an ordinary
        % sum-squared-error calculation)
        chi2 = @(p) sum( ((yavg - fit(p, x)) ./    1).^2 );
        
    otherwise
        error('Chi^2 method not recognized');
end
 
% Initial condition for fminsearch 
x0 = [A b];
 
% Use 'fminsearch' to find the parameters that minimize the chi2 function
pmin = fminsearch(chi2, x0);
 
% Find minimum chi2 using parameters
minchisq = chi2(pmin);
 
% Calculate the "goodness of fit" R^2
Rsq_chi2 = 1 - sum((yavg - fit(pmin, x)).^2)/sum((yavg - mean(yavg)).^2);
 
%% Chi-squared min + 1 
 
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
 
 
%% OLS and WLS via matrix methods
 
% If data is not linear, linearize it!
% NOTE that this will not work for the sine example... 
x_lin = x;
y_lin = yavg; 
 
% OLS
 
% Construct x matrix with one column for each parameter in the model.
% Each column is scaled to the order of the polynomial of x. For a linear
% equation, the first column is x^1 (just x), and second column is x^0,
% which is a vector of ones. This method works for n-degree polynomial fits
 
x_ols = [x_lin, ones(length(x_lin), 1)];
 
% Matrix version: (X'*X)^(-1) * X' * Y
% p_ols = inv(x_ols' * x_ols) * x_ols' * yavg;
 
% But wait! Matlab has a convenient shorthand for this calculation 
% using the backslash '\' operator:
p_ols = x_ols \ y_lin;
 
% OLS coefficients
A_ols = p_ols(1);
b_ols = p_ols(2);
 
% Goodness of fit
Rsq_ols = 1 - sum((y_lin - fit(p_ols, x)).^2)/sum((y_lin - mean(y_lin)).^2);
 
 
% WLS
 
% x and y get weighted by their uncertainty
x_wls = x_ols ./ ysem;
y_wls = y_lin ./ ysem;
 
% Compute coeficients (same as OLS method)
p_wls = x_wls \ y_wls;
 
% Covaraince matrix
cov = inv(x_wls'*x_wls);
 
% NOTE: this this is equivalent to a 1/dy^2 weighting, which can be  
% calculated explicitly by defining a weight matrix of 1/dy^2 on the
% diagonal (zeros everywhere else), then computing a solution using the
% OLS matrices multiplied by the weights like so:
% W = diag(1 ./ (ysem.^2));
 
% % The covariance matrix is
% cov = inv(x_ols' * W * x_ols);
 
% % The coefficients are computed as
% p_wls = cov * x_ols' * W * y_lin;
 
% Try it out yourself and compare. There are a number of other weighting
% strategies that can be employed too, like 1/dy, 1/x, etc... try a few and
% see how they affect the results.
 
 
% Individual coefficients
A_wls = p_wls(1);
b_wls = p_wls(2);
 
% The uncertainties in each parameter lie on the diagonal of the covariance
% matrix
sigma_p = sqrt(diag(cov));
dA_wls = sigma_p(1);
db_wls = sigma_p(2);
 
% Goodness of fit
Rsq_wls = 1 - sum((y_lin - fit(p_wls, x)).^2)/sum((y_lin - mean(y_lin)).^2);
 
 
%% Display
 
% Plot the data model with the minimized parameters
figure('position', [100 100 1500 600]); 
subplot(1,2,1);hold on;
set(gca, 'FontSize', 14);
errorbar(x, yavg, ysem, 'ok');
plot(x, fit(pmin, x), 'k');
 
% Plot OLS/WLS fits if NOT sine example
if ~strcmp(examp, 'sine')
    plot(x, fit(p_ols, x), '--k');
    plot(x, fit(p_wls, x), '-.k');
end
 
xlabel('Domain (Time, Distance, etc..)');
ylabel('Range (Intensity, Voltage, Counts, etc..)');
title(titlestr);
 
% Fit information
str_chi2 = sprintf(['\\chi^2 Fit Params:\n' ...
               'A = %3.2f \\pm %3.2f \n' ...
               'b = %3.2f \\pm %3.2f \n' ...
               'R^2 = %4.3f \n' ...
               '\\chi^2 = %4.3f \n' ...
               '\\chi^2_{red} = %4.3f \n'], pmin(1), dA, pmin(2),...
                db, Rsq_chi2, minchisq, minchisq/dof);
            
str_ols = sprintf(['OLS Fit Params: \n' ...
               'A = %3.2f, b = %3.2f \n' ...
               'R^2 = %4.3f \n'], A_ols, b_ols, Rsq_ols);
           
str_wls = sprintf(['WLS Fit Params: \n' ...
               'A = %3.2f \\pm %3.2f \n' ...
               'b = %3.2f \\pm %3.2f \n' ...
               'R^2 = %4.3f \n'], A_wls, dA_wls, b_ols, db_wls, Rsq_ols);
 
legend({'Data', str_chi2, str_ols , str_wls}, 'Location', 'Best');
 
% Plot the Chi^2 + 1 ellipsoid in parameter space
subplot(1,2,2); hold on;
set(gca, 'FontSize', 14);
contourf(X, Y, Z, 50);
scatter(Arange(idxA), brange(idxb), 'mo', 'MarkerFaceColor', 'Magenta');
plot(pmin(1), pmin(2), 'kx', 'MarkerSize', 15);
xlabel('A'); ylabel('b'); ch = colorbar; ylabel(ch, sprintf('\\chi^2'));
title(sprintf('\\chi^2 Parameter Space'));
 
% Chi^2 + 1 uncertainty boundaries
line([dAmin dAmin], [pmin(2)-r pmin(2)+r],     'Color', 'magenta', 'LineStyle', '--', 'LineWidth', 2);
line([dAmax dAmax], [pmin(2)-r pmin(2)+r],     'Color', 'magenta', 'LineStyle', '--', 'LineWidth', 2);
line([pmin(1)-r pmin(1)+r],     [dbmin dbmin], 'Color', 'magenta', 'LineStyle', '--', 'LineWidth', 2);
line([pmin(1)-r pmin(1)+r],     [dbmax dbmax], 'Color', 'magenta', 'LineStyle', '--', 'LineWidth', 2);
text(pmin(1), pmin(2) - r/2, sprintf('dA = %5.4f', dA), 'Color', 'White', 'FontSize', 14);
text(pmin(1)+r/2, pmin(2), sprintf('dw = %5.4f', db), 'Color', 'White', 'FontSize', 14);
 


