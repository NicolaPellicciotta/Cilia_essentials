function [ ] = plots_Iqtau_colloids2(varargin)
%%Performs plots etc.
%Called with six arguments from Batch2DDM/fit_..., with 1 arg. when
%supplying just the 'fitdir', or none to get the 'fitdir' via GUI.

disp('_ _ _ _');
disp('plots_Iqtau_colloids2');
refit = 'y';
count = 0;

%% Parse i/p: accepts all or no arguments
while (count == 0)
defnargs = 7;   % 6 default arguments

switch nargin
    case defnargs
        Iqtau = varargin{1};
        fits_array = varargin{2};
        tau_q_array = varargin{3}; 
        tau = varargin{4}; 
        q = varargin{5}; 
        fitdir = varargin{6};
        fun = varargin{7};
    otherwise  
        if nargin == 1
            fitdir = varargin{1};
        elseif nargin<1
            ipmess = 'Select the folder containing the .mat file storing the desired fits data';
            fitdir = uigetdir('C:\MicroscopeData',ipmess);
        else
            disp('Error fitdir');
        end
    
        dummy = strsplit(fitdir,filesep);          % Split up dirname according to OS
            fitname = dummy{end};                       % Last entry is comp[utation]name
            ext = 'mat';
            ipfilename = strjoin({fitname,ext},'.');
            path = strjoin({fitdir, ipfilename}, filesep);
        variables = {'Iqtau2','fits_array','tau_q_array','tau','q','fun'};
        file = load(path, variables{:});
            Iqtau = file.Iqtau2;
            fits_array = file.fits_array;
            tau_q_array = file.tau_q_array;
            tau = file.tau;
            q = file.q;
            fun = file.fun;
end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Settings
    params = csvread('Settings.txt',1,0); % px2mum, Th, tempStpt, freqStpt, fun, Stptfoo, fixed

    ipfilename = 'vidpars.mat';                    % vid[eo] par[ameter]s
        dummy = strsplit(fitdir,filesep);          % Split up dirname according to OS
        fummy = dummy(1:end-2);                    % Last entry is comp[utation]name
        userdir = strjoin(fummy, filesep);
    path = strjoin({userdir,ipfilename}, filesep);
    list = ls(path);

    if isempty(list)
        px2mum = params(1);     % If first run with this video
    else
        variable = 'px2mum';    % Has been saved if run before
        file = load(path, variable);
        px2mum = file.px2mum;
    end
    Stptfoo = params(6);   % Start points for fit of D
    fixed = params(7);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

count = 1;
end

%% Plot of I(tau) against tau, interactively scrolling through q
scroll_through_Iqtau_fits(Iqtau,fits_array,tau);
% scroll_through_Iqtau_fits_Split(Iqtau,fits_array,tau);
scroll_through_Iqtau_fits_log(Iqtau,fits_array,tau);

%% Parameter plots
if (fun == 1) || (fun == 2) || (fun == 4) || (fun == 5) || (fun == 6)       % For characteristic time
    wait1 = input('Continue to characteristic time [y/n]? ','s');  
    
    while(strcmp(wait1,'y'))
    [Dmin_q_fit, Dmax_q_fit] = plot_Iqtau_Choice(q, tau_q_array, 1);                        % Plot of tau_q vs q: choose fit range
    [ Dfoo, Dfoo2 ] = plot_trend_fit(Dmin_q_fit,Dmax_q_fit, q, tau_q_array, 1, px2mum, Stptfoo);  % Fit trendlines tau_q vs q. foo forced, foo2 fixed
    plot_Iqtau_Results(px2mum, fitdir, Dmin_q_fit, Dmax_q_fit, q, tau_q_array, 1, Dfoo, Dfoo2); % Plot of tau_q vs q
    plot_Iqtau_Save(px2mum, Dmin_q_fit, Dmax_q_fit, 1, Dfoo, Dfoo2, fitdir)                                  % Save figure
    wait1 = 'n';
    end
end

if (fun == 3) || (fun == 2) || (fun == 4) || (fun == 5) || (fun == 6)       % For oscillation frequency 
    wait2 = input('Continue to oscillation frequency [y/n]? ','s');  
    
    while(strcmp(wait2,'y'))
    [vmin_q_fit, vmax_q_fit] = plot_Iqtau_Choice(q, tau_q_array, 3);                        % Plot of omega vs q: choose fit range
    [ vfoo, vfoo2 ] = plot_trend_fit(vmin_q_fit,vmax_q_fit, q, tau_q_array, 3, px2mum, Stptfoo);  % Fit trendlines omega vs q 
    plot_Iqtau_Results(px2mum, fitdir,vmin_q_fit, vmax_q_fit, q, tau_q_array, 3, vfoo, vfoo2); % Plot of omega vs q
    plot_Iqtau_Save(px2mum, vmin_q_fit, vmax_q_fit, 3, vfoo, vfoo2, fitdir)                                  % Save figure
    wait2 = 'n';
    end
end

%% Continue?  Refit over a different range?
while (strcmp(refit,'y'))
    refit = input('Refit trendline over a different range [y/n]? ','s');     
    disp(' ')
    if(strcmp(refit,'y'))
        switch nargin                                                                      % Recalls itself with appropriate 'nargin'
            case 6
                plots_Iqtau_colloids2(Iqtau, fits_array, tau_q_array, tau, q, fitdir);     % Do plots (but not the fits!) again 
                refit = 'n';
            otherwise  
                plots_Iqtau_colloids2(fitdir);
                refit = 'n';
        end
    end
end

disp('_____________________________');
disp(' ');

end