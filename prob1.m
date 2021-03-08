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


%Part A
%convert to an array to do math
yarray = table2array(rawydata);
%transpose this so that it's in the same format as rawxdata
ymean = mean(yarray, 2);
%we konow the number of samples is Nrep 
yse = (std(yarray')/sqrt(Nrep))';
%let's concatenate this into a table for plotting later
condata = addvars(rawxdata, ymean, yse);

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


figure('Renderer', 'painters', 'Position', [10 10 900 600])
errorbar(condata.x, condata.ymean, condata.yse, 'x');
hold on
plot(x, fit(pmin, x), 'k');
hold off
xlim([0 2.1]);
xticks(0:.1:2.1);
ylim([0 45]);


%temporarily set dA = 1 and db = 1
dA = 1;
db = 1;
str_chi2 = sprintf(['\\chi^2 Fit Params:\n' ...
               'A = %3.2f \\pm %3.2f \n' ...
               'b = %3.2f \\pm %3.2f \n' ...
               'R^2 = %4.3f \n' ...
               '\\chi^2_{min} = %4.3f \n' ...
               '\\chi^2_{red} = %4.3f \n'], pmin(1), dA, pmin(2),...
                db, Rsq_chi2, minchisq, redchisq);
legend('Data', sprintf('\\chi^2_{min}'), 'Location','Northwest');

a = gca; % get the current axis;
% set the width of the axis (the third value in Position) 
% to be 60% of the Figure's width
a.Position(3) = 0.7;            
annotation('textbox', [.85, .5, .1,.1], 'String',  str_chi2); 



%Part B

%Part C

%Part D