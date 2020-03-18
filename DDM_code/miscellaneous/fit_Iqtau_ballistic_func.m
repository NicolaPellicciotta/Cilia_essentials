function [ballistic_fits] = fit_Iqtau_ballistic_func(Iqtau)
%% fit_Iqtau_brownian_func fits the ISF (DDM output) to find A(q), B(q) and \tau_c(q)

%% set general parameters

frame_size = 2*size(Iqtau,1);
max_mode = size(Iqtau,1);      %trying to extract data from the whole Iqtau matrix. Do not change max_mode though or q will be wrong in next script.
max_lag = size(Iqtau,2);


%% set variables, fit functions

% indipendent variable
tau = (1:max_lag)'; % column because of fit

% function to fit
ft = fittype('a*(1-  (sin(tau/tau_c))/(tau/tau_c)  )+b','independent','tau');


%% for loop on the modes, fitting each row if Iqtau with ft

% for mc = max_mode : -1 : 1 %sneaky allocation
parfor mc = 1:max_mode %sneaky allocation
    
    % row of the Iqtau to be fitted
    y = Iqtau(mc,:)'; % column because of fit
    
    % normalising y to make fit work better
    max_y = max(y);
    y = y./ max_y;
    
    %defining fit options
    fo = fitoptions('Method','NonLinearLeastSquares');
    fo.StartPoint = [1, 0, 30];   % options for the fit
    fo.Lower = [0, 0, 0]; %in order: a, b, tau_c
    fo.Upper = [Inf, y(1), Inf];
    
    % actual fit
    ballistic_fits(mc).fit_out = fit(tau,y,ft,fo);
    
    % save normalisation factor so we can use it for plots later
    ballistic_fits(mc).norm_fact = max_y; % We'll multiply a and c by it, or divide Iqtau(mc,:)
    
end %for

end