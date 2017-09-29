%% XY to Hough Space Demonstration
% David Olson
% 26 Sep 17

close all;
clear all;
clc;

%% Define Points in Space
x = 2 : 0.5 : 7;
y = 2 : 0.5 : 7;

xyPoints = [x; y];

%% Define Theta Bins and Empty Hough 2-D Histogram

theta = (-pi/2) : (0.5 * (pi/180)) : (pi/2);

% Empty Hough Histogram Size
M = length(-10 : 0.1 : 10);      % # of Rows
N = length(theta);                % # of Columns

% Allocate Memory for 2-D Histogram
HoughHist = zeros(M, N);  % M x N matrix

%% Begin Iteration of Theta with Given Points

P_1__t = [0; 15];
P_1__b = [0; -15];


for k = 1 : length(xyPoints)
    
    for n = 1 : N
        
        rho = xyPoints(1,k)*cos(theta(n)) + xyPoints(2,k)*sin(theta(n));
        rowIndex = ((round(rho, -log10(0.1)) * (1/0.1)) + 1) + (10*(1/0.1));
        HoughHist(rowIndex, n) = HoughHist(rowIndex, n) + 1;
        C_0__1 = [cos(theta(n)), -sin(theta(n)); sin(theta(n)) cos(theta(n))];
        P_0__t = (C_0__1 * P_1__t) + [x(k); y(k)];
        P_0__b = (C_0__1 * P_1__b) + [x(k); y(k)];
        
        % For First Iteration
        if ((k == 1) && (n == 1))
            subplot(1, 2, 1)
            plot(xyPoints(1,:), xyPoints(2,:), 'bs')
            hline = line([P_0__b(1), P_0__t(1)], [P_0__b(2), P_0__t(2)], 'LineWidth', 2, 'Color', 'r');
            rholine = line([0, rho*cos(theta(n))], [0, rho*sin(theta(n))], 'LineWidth', 1, 'Color', 'g');
            line([0, 0], [-5 10], 'LineWidth', 0.5, 'Color', 'k')
            line([-5, 10], [0 0], 'LineWidth', 0.5, 'Color', 'k')
            title('X-Y Space')
            xlabel('X axis')
            ylabel('Y axis')
            xlim([-5 10])
            ylim([-5 10])
            grid on
        else
            set(hline, 'XData', [P_0__b(1), P_0__t(1)]);
            set(hline, 'YData', [P_0__b(2), P_0__t(2)]);
            set(rholine, 'XData', [0, rho*cos(theta(n))]);
            set(rholine, 'YData', [0, rho*sin(theta(n))]);
        end
        
        % Plot 2-D Histogram
        if ((k == 1) && (n == 1))
            subplot(1, 2, 2)
            hmesh = mesh(HoughHist);
            title('Hough Space')
            xlabel('Theta')
            ylabel('Rho')
            view(0, 90)
            xlim([0 N])
            ylim([0 M])
        else
            set(hmesh, 'CData', HoughHist);
        end
        
        pause(0.001)
        
    end
end