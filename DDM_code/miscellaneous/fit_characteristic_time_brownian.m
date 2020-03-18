function[fit_out_brownian, fit_out_freexp] = fit_characteristic_time_brownian(brownian_fits)
%% fit_characteristic_time_brownian fits tau_c(q) to find the diffusion coefficient

%% set parameters

max_mode = numel(brownian_fits);

% change the next two lines if you want to fit a different range (in q) of the datapoints
min_mode_fit = 1;
max_mode_fit = max_mode;


%% set variables, fit functions


% indipendent variable
q = pi/max_mode*(1:max_mode_fit)'; % column because of fit

% functions to fit
ft_brownian = fittype('1/(Dm*q^2)','independent','q');
ft_freeexp  = fittype('1/(Dm*q^alpha)','independent','q');


%% extract data from brownian_fits structure

modes       = (1:max_mode)';     % column because of fit
tau_c       = zeros(max_mode,1); % column because of fit
err_tau_c   = zeros(max_mode,1); % column because of fit

% using a for loop is more readable than using arrayfun
for i = modes'
    
    % characteristic time
    tau_c(i) = brownian_fits(i).fit_out.tau_c;
    
    % the fit's confidence interval
    confs = 0.5*diff(confint(brownian_fits(i).fit_out,0.682));
    err_tau_c(i) = confs(3);
    
end %for


%% call function to set interval of modes to fit

[hf, min_mode_fit, max_mode_fit] = fit_characteristic_time_choose_range(modes, tau_c, err_tau_c);
waitfor(hf)
min_mode_fit, max_mode_fit
%select range on arrays
q_tofit         =         q(min_mode_fit:max_mode_fit);
tau_c_tofit     =     tau_c(min_mode_fit:max_mode_fit);
err_tau_c_tofit = err_tau_c(min_mode_fit:max_mode_fit);


%% fit of the time constant tau_c vs q

hf = figure;

ha = axes;
hold on;
box on;
ha.XScale = 'log';
ha.YScale = 'log';

% plot all the curve
he1 = errorbar(q, tau_c, err_tau_c);
he1.Marker = '.';
he1.MarkerSize = 8;
he1.Color = 'b';
he1.LineStyle = 'none';
he1.LineWidth = 1.2;

% plot only the selected interval
he2 = errorbar(q_tofit, tau_c_tofit, err_tau_c_tofit);
he2.Marker = '.';
he2.MarkerSize = 8;
he2.Color = 'g';
he2.LineStyle = 'none';
he2.LineWidth = 1.2;

%set labels
hxl = xlabel('Spatial Frequency q, [px^{-1}]', 'FontSize',16);
hyl = ylabel('Characteristic time \tau_c, [frames]', 'FontSize',16);


%% actual fit

% fit options
fo_brownian = fitoptions('Method','NonLinearLeastSquares','Weights',1./err_tau_c_tofit.^2,'StartPoint',[1]);
fo_freexp = fitoptions('Method','NonLinearLeastSquares','Weights',1./err_tau_c_tofit.^2,'StartPoint',[1 2]);

% fit with the two equations
fit_out_brownian = fit(q_tofit,tau_c_tofit,ft_brownian,fo_brownian);
fit_out_freexp   = fit(q_tofit,tau_c_tofit,ft_freeexp ,fo_freexp);


%% plot fit on graph

% plot stuff
hpf = plot(q_tofit,fit_out_freexp(q_tofit)  ,'k', 'LineWidth', 1.2);

hpb = plot(q_tofit,fit_out_brownian(q_tofit),'r', 'LineWidth', 1.2);

ha.XScale = 'log';
ha.YScale = 'log';

hleg = legend('non-fitted data points',...
    'fitted data points',...
    ['fit with free exponent \alpha = ',sprintf('%.2f',fit_out_freexp.alpha)],...
    'fit with exponent fixed at 2',...
    'Location','SouthWest');



end

