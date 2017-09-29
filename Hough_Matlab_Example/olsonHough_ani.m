function [m_meas, c_meas, HoughHist] = olsonHough(xyPoints, thetaInc, rhoRange, rhoInc)
%OLSONHOUGH This function is a test of a self-made Hough Transform

% xyPoints -- 2xN array containing x-y points
% thetaInc -- Range is predefined -90 to 90 deg, increment needs input in
% rad
% rhoRange -- range of rho to search for, needs positive number
% rhoInc -- increment of rho, must be a multiple of ten!  i.e. 0.1, 0.01,
% etc.

%% Define Theta Bins and Empty Hough 2-D Histogram

theta = (-pi/2) : thetaInc : (pi/2);

% Empty Hough Histogram Size
M = length(-rhoRange : rhoInc : rhoRange);      % # of Rows
N = length(theta);                              % # of Columns

% Allocate Memory for 2-D Histogram
HoughHist = zeros(M, N);  % M x N matrix

%% Gather Votes for Possible Lines within the x-y provided points

for k = 1 : length(xyPoints)
    
    for n = 1 : N
        
        rho = xyPoints(1,k)*cos(theta(n)) + xyPoints(2,k)*sin(theta(n));
        rowIndex = ((round(rho, -log10(rhoInc)) * (1/rhoInc)) + 1) + (rhoRange*(1/rhoInc));
        HoughHist(rowIndex, n) = HoughHist(rowIndex, n) + 1;
        
    end
    
    % Plot XY Data Points
    subplot(1, 2, 1)
    hold on
    plot(xyPoints(1,:), xyPoints(2,:), 'bs')
    hold on
    plot(xyPoints(1,k), xyPoints(2,k), 'r*')
    hold on
    plot(xyPoints(1,(1:k-1)), xyPoints(2,(1:k-1)), 'y*')
    title('XY Points Evaluation')
    xlabel('X axis')
    ylabel('Y axis')
    xlim([-5 5])
    ylim([-5 10])
    grid on
    
    % Plot 2-D Histogram
    subplot(1, 2, 2)
    mesh(HoughHist)
    title('2-D Histogram of Rho and Theta Parameters')
    xlabel('Theta Bins')
    ylabel('Rho Bins')
    view(0, 90)
    xlim([0 361])
    ylim([0 401])
    
    pause(0.5)
    
end

%% Determine Peak of 2-D Histogram and best choice for Rho and Theta

% Determining Bin Locations of the peak of the 2-D Histogram
[maxOfThetas, thetaMaxes] = max(HoughHist);
[highestVote, bestBinTheta] = max(maxOfThetas);
bestBinRho = thetaMaxes(bestBinTheta);

% Determine Calculated Rho's and Theta's
rho_meas = -rhoRange + (bestBinRho * rhoInc);
theta_meas = (-pi/2) + (bestBinTheta * thetaInc);

%% Determine Measured 'm' and 'c'

m_meas = -cot(theta_meas);
c_meas = rho_meas * (sin(theta_meas) - m_meas*cos(theta_meas));

return

end

