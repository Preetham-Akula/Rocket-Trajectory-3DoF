function [fitresult, gof] = rho_h_CurveFit(h, rho)%,fig_num)
%CREATEFIT(H,RHO)
%  Create a fit.
%
%  Data for 'Rho_h' fit:
%      X Input : h
%      Y Output: rho
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 02-Feb-2020 11:41:21


%% Fit: 'Rho_h'.
[xData, yData] = prepareCurveData( h, rho );

% Set up fittype and options.
ft = fittype( 'fourier5' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0 0 0 0 0 0 0 0 0 0 0 0.0261799387799096];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
% figure( 'Name', 'Rho_h' );
figure;%(fig_num)
hold on;
h_1 = plot( fitresult, xData, yData );
legend( h_1, 'rho vs. h', 'Rho_h', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'h', 'Interpreter', 'none' );
ylabel( 'rho', 'Interpreter', 'none' );
grid on


