function [ ] = fit_Iqtau_colloids2(varargin)
% fit_Iqtau_colloids.m fits exponentials to find the characteristic times, tau_q, and then 'calls' 
% refit_Iqtau_colloids.m to produce the choice scroll, choice, and end plots fun frames2s, Iqtau, DDMdir

%%
disp('_ _ _ _ ');
disp('fit_Iqtau_collloids2');

%% Parse i/p: accepts all or no arguments
defnargs = 4;       % def[ault] n[umber of] arg[ument]s

switch nargin
    case defnargs                   % Run from BatchDDM2
        fun = varargin{1};
        frames2s = varargin{2};
        Iqtau = varargin{3};
        DDMdir = varargin{4};        
    otherwise                       % Run from command line
        if nargin == 1              
            DDMdir = varargin{1};
        elseif nargin<1
            ipmess = 'Select the folder containing the .mat file storing the desired Iqtau matrix';
            DDMdir = uigetdir('C:\MicroscopeData',ipmess);
        else
            disp('Error DDMdir');
        end
        
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Settings
            params = csvread('Settings.txt',1,0); % px2mum, Th, tempStpt, freqStpt, fun, Stptfoo, fixed

            tempStPt = params(3); % Automatic value should work.  Run again without
            freqStPt = params(4);    % arguments to try alternative, if necessary
            fun = params(5);            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        dummy = strsplit(DDMdir,filesep);          % Split up dirname according to OS
        ext = 'mat';
            fummy = dummy(1:end-1);                       % Last entry is comp[utation]name
            userdir = strjoin(fummy, filesep);
            ipfilename = strjoin({'vidpars',ext},'.');
        path = strjoin({userdir, ipfilename}, filesep);
        variables = {'frames2s'};
        file = load(path, variables{:});
        frames2s = file.frames2s;

            compname = dummy{end};                       % Last entry is comp[utation]name
            ipfilename = strjoin({compname,ext},'.');
        path = strjoin({DDMdir, ipfilename}, filesep);
        variables = {'Iqtau'};
        file = load(path, variables{:});
        Iqtau = file.Iqtau;
end

tic;        % Start timing

%% Parameters
Iqtau2 = 1e4 .* Iqtau;          % To make the fit work, c.f. discretisation of numbers
frame_size = 2*size(Iqtau2,1);
max_q = size(Iqtau2,1);         % Attempt to extract data from the whole Iqtau matrix 
max_tau = size(Iqtau2,2);       % Tighten these 2 parameters if the fit gives errors
frames = (1:max_tau)';
tau = frames2s*frames;          % lag times in seconds in correct dimension for later

% Start points for fits of parameters of saturation functions
if nargin<defnargs                  % If called from CL, SPs taken from 'Settings.txt'
    disp(['The set start point for the characteristic time is ', num2str(tempStPt),'.']);
    disp(['The set start point for the oscillation frequency is ', num2str(freqStPt),'.']);
else
    tempStPt = tau(end)/10;         % If run from BatchDDM2, SPs generated here
    disp(['The automatic start point for the characteristic time is ', num2str(tempStPt),'.']);
    if(fun~=1)
        freqStPt = 10;
        disp(['The autoatic start point for the oscillation frequency is ', num2str(freqStPt),'.']);
    end
end

%% Fit
% Which function?
switch fun          % Define the functions that can be fit, func for fittype, funname for titles
    case 1 
        func = 'a*(1-exp(-tau/b))+c';
        funname = 'exp';
    case 2
        func = 'a*(1-exp(-tau/b)*sinc(d*tau))+c';
        funname = 'expsinc';
    case 3
        func = 'a*(1-sinc(d*tau))+c';
        funname = 'sinc';
    case 4
        func = 'a*(1-exp(-(tau/b)^2)*sinc(d*tau))+c';
        funname = 'gausssinc';
    case 5
        func = 'a*(1-exp(-tau/b)*cos(d*tau))+c';
        funname = 'expcos';
    case 6
        func = 'a*(1-exp(-(tau/b)^2)*cos(d*tau))+c';
        funname = 'gausscos';
end

% Initialise arrays
qmode = (1:max_q)';
q = 2*pi/frame_size*qmode; %in px^-1 % definition of the array of q
fits_array = cell(max_q, 1);  % Preallocation not helpful/correct for cell array? % A = zeros(max_q,1);   % initialisation of variables for the fit

if (fun == 1)
    tau_q_array = zeros(max_q, 2);  % tau_q = zeros(max_q,1); % B = zeros(max_q,1);% err_a = A;% err_tau_q = tau_q; % err_B = B;
else
    tau_q_array = zeros(max_q, 4);  
end

% Fit
disp(['The function to be fit, number ',num2str(fun),', is f(tau)=', func]); % Function to fit 
ft = fittype(func,'independent','tau'); 

for q_c=1:max_q     % q_c == counter on q
    
    y =  Iqtau2(q_c,1:max_tau)';        % Here extracting a row (function of tau) of the Iqtau matrix     
                                        % Nearly the same fit with three different o/p depending on functional form
    if (fun == 1)                               % With exp. time cst. only                      
            fo = fitoptions('Method','NonLinearLeastSquares','Lower',[0, 0, 0], ...
                    'Upper',[Inf, Inf, y(1)],'StartPoint', ...
                    [max(Iqtau2(q_c,:)'), tempStPt, 0]); %,'Weights',1./err_y.^2);      % options for the fit
            fit_out = fit(tau,y,ft,fo);                     % A(q_c) = fit_out.a;          % saturation value
            tau_q_array(q_c, 1) = fit_out.b;                % tau_q(q_c) = fit_out.b;      % time constant (this is what this fit is done for)
            temp = 0.5*diff(confint(fit_out,0.682));        % B(q_c) = fit_out.c;          % offset
            tau_q_array(q_c, 2) = temp(2);                  % temp and next three lines: errors on A, tau_q and B         
                                                            % err_A(q_c) = temp(1); % err_tau_q(q_c) = temp(2); err_B(q_c) = temp(3);     
    elseif  (fun == 2) || (fun == 4) || (fun == 5) || (fun == 6)        % With exp. time cst. and oscillation frequency
            fo = fitoptions('Method','NonLinearLeastSquares','Lower',[0, 0, 0, 0], ...
                    'Upper',[Inf, Inf, Inf, y(1)],'StartPoint', ...
                    [max(Iqtau2(q_c,:)'), tempStPt, tempStPt, 0]); %,'Weights',1./err_y.^2);      % options for the fit
            fit_out = fit(tau,y,ft,fo);                
            tau_q_array(q_c, 1) = fit_out.b;            % time constant (this is what this fit is done for)
            tau_q_array(q_c, 3) = fit_out.d;            % sinc frequency
            temp = 0.5*diff(confint(fit_out,0.682));     
            tau_q_array(q_c, 2) = temp(2);
            tau_q_array(q_c, 4) = temp(3);  
    elseif (fun == 3)                           % With oscillation frequency only
            fo = fitoptions('Method','NonLinearLeastSquares','Lower',[0, 0, 0], ...
                    'Upper',[Inf, Inf, y(1)],'StartPoint', ...
                    [max(Iqtau2(q_c,:)'), freqStPt, 0]); %,'Weights',1./err_y.^2);      % options for the fit
            fit_out = fit(tau,y,ft,fo);                     
            tau_q_array(q_c, 3) = fit_out.d;                    % sinc frequency (always start at index = 3)
            temp = 0.5*diff(confint(fit_out,0.682));      
            tau_q_array(q_c, 4) = temp(2);         
    end
    fits_array{q_c} = fit_out;
end                

toc         % Stop timing, print
                                                    
%% Save Fit variables
% Save under userdir/DDM_core2_[Th]_[mmdd_HHMM]/fit_[max_q]_[max_tau]_[func]_[mmdd_HHMM]/fit_[max_q]_[max_tau]_[func].mat
t = datestr(now, 'mmdd_HHMM');

fitname = ['fit_',funname,'_',num2str(max_q),'_',num2str(max_tau),'_',num2str(tempStPt),'_',t];
fitdir = strjoin({DDMdir,fitname}, filesep);
if ~isdir(fitdir)
    mkdir(fitdir);
end

ext = 'mat';
filename = strjoin({fitname,ext},'.');
path = strjoin({fitdir,filename}, filesep);
variables = {'Iqtau2','fits_array','tau_q_array','tau','q','fun'};  %,'tau_q_array','q','tau'};

save(path,variables{:});

%% Call the plots function
plots_Iqtau_colloids2(Iqtau2, fits_array, tau_q_array, tau, q, fitdir, fun);

end