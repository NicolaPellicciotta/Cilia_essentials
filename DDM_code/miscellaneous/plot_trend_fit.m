function [ foo, foo2 ] = plot_trend_fit(min_q_fit,max_q_fit, q, tau_q_array, index, px2mum, Stptfoo)
%% Fits trendline within specified range

% Which fit?
switch index
    case 1
        ft = '1/(A*x^2)';
        ft2 = '1/(A*x^alpha)';
        Conv = px2mum^2;
        What = 'Diffusion Coefficient';
        Unit = 'um2/s';
    case 3
        ft = 'A*x';
        ft2 = 'A*x^alpha';
        Conv = px2mum;
        What = 'Velocity';
        Unit = 'um/s';
end

% Do fit
foo = fit(q(min_q_fit:max_q_fit), ...
            tau_q_array(min_q_fit:max_q_fit,index), ...
            fittype(ft,'independent','x'), ...
            fitoptions('Method','NonLinearLeastSquares','StartPoint',Stptfoo, ...
            'Weights',1./tau_q_array(min_q_fit:max_q_fit,(index+1)).^2));
foo2 = fit(q(min_q_fit:max_q_fit),...
            tau_q_array(min_q_fit:max_q_fit,index), ....
            fittype(ft2,'independent','x'),...
            fitoptions('Method','NonLinearLeastSquares','StartPoint',[Stptfoo, 2], ...
            'Weights',1./tau_q_array(min_q_fit:max_q_fit,(index+1)).^2));

% Conversion and print
disp([What,' (Fitted) = ', num2str(foo2.A*Conv), ' ', Unit]);  % now in \mum^2/s (foo.A is in px^2/s)
disp(['Fitted exponent = ', num2str(foo2.alpha),'. ']);
disp([What,' (Forced) = ', num2str(foo.A*Conv),' ', Unit]);
disp(' ');

end 