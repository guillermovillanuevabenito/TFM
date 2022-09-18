function [fitresult, gof] = Network_fit(t,y)
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

ft = fittype( 'rate_model_network_fit(x,tau2e,tau2i,tau3i,x2ethr,x2eslp,x2ithr,x2islp,x3ithr,x3islp)', 'coefficients', {'tau2e','tau2i','tau3i','x2ethr','x2eslp','x2ithr','x2islp','x3ithr','x3islp'},'independent', 'x','dependent', 'y');

opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [0 0 0 -0.5 0 -0.2 0 -0.2 0];
%opts.StartPoint = [400 200 100 -0.2 1.4 0 1.4 0 0.3];
opts.StartPoint = [600 200 100 -0.2 1.4 0 1.4 0 0.3];
opts.Upper = [1000 1000 1000 0 4 0.2 4 0.1 0.5];

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