%% olsonHoughTest m-file
% David Olson
% 6 Sep 2017

close all;
clear all;
clc;

%% Load Data and Provide Function Inputs

load('houghTestData.mat') % lets work with the same data each time

xyPoints = [x_noise; y_noise];      % Provide array of xy points
thetaInc = 0.5 * (pi/180);          % Arguement needs to be in rad
rhoRange = 20;                      % Depends on size of line to find
rhoInc = 0.1;                       % Needs to be multiple of 10 

%% Plot Inital Data

figure
hold on
htrue = plot(x_true, y_true, 'r');
set(htrue, 'LineWidth', 3);
plot(x_noise, y_noise, 'bs')
title('Starting Data:  Truth and Noisy Measurements')
xlabel('X Axis')
ylabel('Y Axis')
grid on
legend('Truth', 'Measurements', 'Location', 'SouthEast')

disp('Hit spacebar to continue')
pause

%% Call Function

[m_meas, c_meas, HoughHist] = olsonHough_ani(xyPoints, thetaInc, rhoRange, rhoInc);

%% Plot Results

x_meas = x_true;
y_meas = (m_meas*x_meas) + c_meas;

m_error = ((m_meas - m_true) / m_true) * 100;
c_error = ((c_meas - c_true) / c_true) * 100;

% Plot mesh
subplot(1, 2, 2)
mesh(HoughHist)
title('2-D Histogram of Rho and Theta Parameters')
xlabel('Theta Bins')
ylabel('Rho Bins')
view(0, 90)
xlim([0 361])
ylim([0 401])

% Plot hough estimate
subplot(1, 2, 1)
htrue = plot(x_true, y_true, 'r');
set(htrue, 'LineWidth', 3);
hnoise = plot(x_noise, y_noise, 'bs');
hmeas = plot(x_meas, y_meas, 'g');
set(hmeas, 'LineWidth', 3);
%legend('Truth', 'Measurements', 'Hough Estimate', 'Location', 'SouthEast')
title('Truth and Measurements and Resultant Line')
xlabel('X axis')
ylabel('Y axis')
xlim([-5 5])
ylim([-5 10])
text(-4, 9.5, ['Measured Slope:  ', num2str(m_meas)]);
text(-4, 9.0, ['Measured Y-Intercept:  ', num2str(c_meas)]);
text(-4, 8.5, ['Slope Percent Error:  ', num2str(m_error), '%']);
text(-4, 8.0, ['Y-Intercept Percent Error:  ', num2str(c_error), '%']);
grid on

