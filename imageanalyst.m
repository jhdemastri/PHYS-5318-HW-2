% Uses fitnlm() to fit a non-linear model (an power law curve) through noisy data.
% Requires the Statistics and Machine Learning Toolbox, which is where fitnlm() is contained.
% Initialization steps.
clc;    % Clear the command window.
close all;  % Close all figures (except those of imtool.)
clear;  % Erase all existing variables. Or clearvars if you want.
workspace;  % Make sure the workspace panel is showing.
format long g;
format compact;
fontSize = 20;
X = [30,29,21,15,9,1];
Y = [0.05,2,3.4,6.18,7.56,10];
% Now we have noisy training data that we can send to fitnlm().
% Plot the noisy initial data.
plot(X, Y, 'b*', 'LineWidth', 2, 'MarkerSize', 15);
grid on;
% Convert X and Y into a table, which is the form fitnlm() likes the input data to be in.
tbl = table(X', Y');
% Define the model as Y = a * (x .^ b)
% Note how this "x" of modelfun is related to big X and big Y.
% x((:, 1) is actually X and x(:, 2) is actually Y - the first and second columns of the table.
modelfun = @(b,x) b(1) * x(:, 1) .^ + b(2);
b = 10 % Close guesses based on what I see from the raw data plot.
m = -1
beta0 = [b, m]; % Guess values to start with.  Just make your best guess.
% Now the next line is where the actual model computation is done.
mdl = fitnlm(tbl, modelfun, beta0);
% Now the model creation is done and the coefficients have been determined.
% YAY!!!!
% Extract the coefficient values from the the model object.
% The actual coefficients are in the "Estimate" column of the "Coefficients" table that's part of the mode.
coefficients = mdl.Coefficients{:, 'Estimate'}
% Now we have the coefficients and we can plot y for ANY x, not just the training set.
% Let's make up a bunch of x (50) from the min to the max.
xFitted = linspace(min(X), max(X), 50);
% Create smoothed/regressed data using the model:
yFitted = coefficients(1) * xFitted .^ coefficients(2);
% Now we're done and we can plot the smooth model as a red line going through the noisy blue markers.
hold on;
plot(xFitted, yFitted, 'r*-', 'LineWidth', 2);
grid on;
title('Power Law Regression with fitnlm()', 'FontSize', fontSize);
xlabel('X', 'FontSize', fontSize);
ylabel('Y', 'FontSize', fontSize);
legendHandle = legend('Noisy Y', 'Fitted Y', 'Location', 'north');
legendHandle.FontSize = 25;
message = sprintf('Coefficients for Y = b * X ^ m:\n  b = %8.5f\n  m = %8.5f',...
  coefficients(1), coefficients(2));
text(5, 3, message, 'FontSize', 23, 'Color', 'r', 'FontWeight', 'bold', 'Interpreter', 'none');
% Set up figure properties:
% Enlarge figure to full screen.
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
% Get rid of tool bar and pulldown menus that are along top of figure.
% set(gcf, 'Toolbar', 'none', 'Menu', 'none');
% Give a name to the title bar.
set(gcf, 'Name', 'Demo by ImageAnalyst', 'NumberTitle', 'Off')