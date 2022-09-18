function [fitresult, gof] = f_alphafit2(t, y)
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
[xData, yData] = prepareCurveData( t, y);

%ft = fittype( 'max(0,a * (x-c)/b .* exp(1-(x-c)/b)) + max(0,a * (x-c-180)/b .* exp(1-(x-c-180)/b)) ', 'independent', 'x', 'dependent', 'y' );

ft = fittype( 'conv_with_kernel(x,a,b,c)', 'independent', 'x','dependent', 'y');

opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [0 1 -100];
opts.StartPoint = [1 400 0];
opts.Upper = [100 10000 1000];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

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

