function [ ] = plot_Iqtau_Results(px2mum, fitdir, min_q_fit, max_q_fit, q, tau_q_array, index, foo, foo2)
%%

%%
% Parse title from directory of fit data
gummy = strsplit(fitdir,filesep); % Split string by file separator ('\' or '/') % Recombine with '_' (for clarity) 
ti = strjoin(strsplit(strjoin(gummy(end-2:end),'_'),'_'),'-'); % Replace by '-' ('_' cause subscripts)

% figure
figure(index);
clf;
hold on;
errorbar(q([1:min_q_fit,max_q_fit:end]), ...
            tau_q_array([1:min_q_fit,max_q_fit:end],index), ...
            tau_q_array([1:min_q_fit,max_q_fit:end],(index+1)), ...
                '.', 'Color', [1,0.64,0]);                      % Plot error bars outside fit range in purple
errorbar(q(min_q_fit:max_q_fit), ...
            tau_q_array(min_q_fit:max_q_fit,index), ...
            tau_q_array(min_q_fit:max_q_fit,(index+1)), ...
                'g.');                                          % Plot error bars inside fit range in green
plot(q(1:end),tau_q_array(1:end,index),'+','MarkerSize',1);     % Plot points in blue

switch index
    case 1
        % plot trendline fits
        plot(q(min_q_fit:max_q_fit),1./(foo.A*q(min_q_fit:max_q_fit).^2),'r');           % fixed (=2) exponent
        plot(q(min_q_fit:max_q_fit),1./(foo2.A*q(min_q_fit:max_q_fit).^foo2.alpha),'k'); % fitted exponent

        % Dressing
        Sca = 'log';
        yla = '\tau_q / s';
        What = 'D';
        Conv = px2mum^2; 
        unit = '\mum^2/s';
        fixed = 2;
        Loc = 'SouthWest';         
    case 3
        % plot trendline fits
        plot(q(min_q_fit:max_q_fit),(foo.A*q(min_q_fit:max_q_fit)),'r');           % fixed (=1) exponent
        plot(q(min_q_fit:max_q_fit),(foo2.A*q(min_q_fit:max_q_fit).^foo2.alpha),'k'); % fitted exponent

        % Dressing
        Sca = 'log';
        yla = '\omega / rads^{-1}';
        What = 'v';
        Conv = px2mum;
        unit = '\mum/s';
        fixed = 1;
        Loc = 'NorthWest';
end

% Dressing
set(gca,'XScale',Sca,'YScale',Sca); % log-log plot
xlabel('q / px^{-1}','FontSize',14);    % q now in (px)^-1 (*not* mode number)
ylabel(yla,'FontSize',14);     % characteristic time in s   
title(ti, 'FontSize', 14);              % Title from above
legend('non-fitted data points'' errorbars','fitted data points'' errorbars', ...
        ['q: ', num2str(min_q_fit),'-',num2str(max_q_fit),' ', What,' = ',num2str(foo2.A*Conv),' ', unit],...
        ['fit with exponent fixed at ', num2str(fixed)], ['fit with free exponent \alpha = ',num2str(foo2.alpha)], ...
        'Location',Loc);        % Legend displays fit range, trendline exponents, and velocity 
end
