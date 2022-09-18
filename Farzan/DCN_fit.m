function [fitresult, gof] = DCN_fit(t, y)
%CREATEFIT(ILIM,ERRL,ILIN)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : t
%      Y Output: y
%      Weights : weights
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%% Fit:
[xData, yData] = prepareCurveData(t, y);

ft = fittype( 'rate_fit(x,tau1,tau2,m1,k1,m2,k2)', 'coefficients', {'tau1','tau2','m1','k1','m2','k2'},'independent', 'x','dependent', 'y');

opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [0 0 0 0 0 0];
opts.StartPoint = [1000 500 0.3 0.15 0.1 0.2];
opts.Upper = [2000 1000 1 1 1 1];

% Fit model to data.
[fitresult, gof] = fit(xData, yData, ft, opts);

%%
% Plot fit with data.
% figure( 'Name', 'untitled fit 1' );
% h = plot( fitresult, xData, yData );
% legend( h, 'errL vs. ilim with iliN', 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
% % Label axes
% xlabel( 'ilim', 'Interpreter', 'none' );
% ylabel( 'errL', 'Interpreter', 'none' );
% grid on
% 