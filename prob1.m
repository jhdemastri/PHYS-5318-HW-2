%Author: John DeMastri
%Date: 3/9/2021
%For: PHYS 5318 HW #2

clear; close all;
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
pts  = 14;       % Data points 
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
% If data is not linear, linearize it!
% NOTE that this will not work for the sine example... 
x_lin = log(x);
y_lin = log(ymean); 
 
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
x_wls = x_ols ./ yse;
y_wls = y_lin ./ yse;
 
% Compute coeficients (same as OLS method)
p_wls = x_wls \ y_wls;
 
% Covaraince matrix
cov = inv(x_wls'*x_wls);

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

%% Plot and Output
figure('Renderer', 'painters', 'Position', [10 10 900 600])
errorbar(condata.x, condata.ymean, condata.yse, 'x');
hold on
plot(x, fit(pmin, x), 'k');
plot(x, fit(p_ols, x), '--k');
plot(x, fit(p_wls, x), '-.k');
hold off
xlim([0 2.1]);
xticks(0:.1:2.1);
ylim([0 45]);


%temporarily set dA = 1 and db = 1
%dA = 1;
%db = 1;
str_chi2 = sprintf(['\\chi^2 Fit Params:\n' ...
               'A = %3.2f \\pm %3.2f \n' ...
               'p = %3.2f \\pm %3.2f \n' ...
               'R^2 = %4.3f \n' ...
               '\\chi^2_{min} = %4.3f \n' ...
               '\\chi^2_{red} = %4.3f \n'], pmin(1), dA, pmin(2),...
                db, Rsq_chi2, minchisq, redchisq);
            
str_ols = sprintf(['OLS Fit Params: \n' ...
               'A = %3.2f, p = %3.2f \n' ...
               'R^2 = %4.3f \n'], A_ols, b_ols, Rsq_ols);
           
str_wls = sprintf(['WLS Fit Params: \n' ...
               'A = %3.2f \\pm %3.2f \n' ...
               'p = %3.2f \\pm %3.2f \n' ...
               'R^2 = %4.3f \n'], A_wls, dA_wls, b_ols, db_wls, Rsq_ols);
legend('Data', sprintf('\\chi^2_{min}'), 'ols', 'wls', 'Location','Northwest');

a = gca; % get the current axis;
% set the width of the axis (the third value in Position) 
% to be 60% of the Figure's width
a.Position(3) = 0.7;            
annotation('textbox', [.84, .825, .1,.1], 'String',  str_chi2); 
annotation('textbox', [.84, .525, .1,.1], 'String',  str_ols);
annotation('textbox', [.84, .225, .1,.1], 'String',  str_wls);