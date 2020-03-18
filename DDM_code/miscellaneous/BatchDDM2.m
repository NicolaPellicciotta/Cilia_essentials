function [ ] = BatchDDM2(varargin)
%%Runs sequence of functions for DDM analysis (v2)
%Usage: BatchDDM2() accepts up to three arguments: 
%-- 's' accept settings in settings.txt
%-- '[string == ['C:\*'] | ['/*']]' directory of video
%-- '[number]' framerate

%% Check I/p
% Initialise
check = 'x';
framerate = 0;  % Can't have zero framerate c.f. later test
userdir = '';   %  ''    ''  empty directory c.f.   ''
settings = '';  % empty var c.f. later test

% Read in inputs given
for i = 1.0:1:nargin            
    if(isfloat(varargin{i}))
        framerate = varargin{i};        % Only non-default float is framrate
        frames2s = framerate;
    elseif(strcmp(varargin{i},'s'))
        settings = 's';                 % Have 's' to use all settings
    elseif(logical(regexp(varargin{i},'\w:\\*')) || logical(regexp(varargin{i},'/\*', 'Once')))  % N.B. \w == alphanumeric character
        userdir = varargin{i};                                         % Directory in Win or Linux
    else 
        disp(['error in argument ', num2str(i),'.'])
    end
end

% If inputs not given
if(isempty(userdir))            % userdir
    dirmess = 'Please select the folder that contains the desired microscope video';
    userdir = uigetdir('C:\MicroscopeData',dirmess);
end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Settings
    params = csvread('Settings.txt',1,0); % px2mum, Th, tempStpt, freqStpt, fun, Stptfoo, fixed
        ipfilename = 'vidpars.mat';
        path = strjoin({userdir,ipfilename}, filesep);
        list = ls(path);
    if isempty(list)
        px2mum = params(1);     % If first run with this video
    else
        variable = 'px2mum';    % Has been saved if run before
        file = load(path, variable);
        px2mum = file.px2mum;
    end
    Th = params(2);
    tempStPt = params(3);
    freqStPt = params(4);
    fun = params(5);
    Stptfoo = params(6);
    fixed = params(7);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(framerate==0)            % frames2s
    if isempty(list)
        framerate = input('Enter the framerate of the video: ');
        frames2s = 1/framerate;
    else
        variable = 'frames2s';
        file = load(path, variable);
        frames2s = file.frames2s;
    end
end

% If flag 's' not given in function call, get 'Settings' here:
if(isempty(settings))
    
    pmess = ['Enter the microscope''s micron/px ratio \n[enter for default=',num2str(px2mum),']:'];
    pm = input(pmess);
    if ~(isempty(pm))
           px2mum = pm;               % px2mum: multiply measure in px by px2mum to have the measure in mum
    end    
    
    hTmess = ['Enter the threshold for throwing away \nframes in the analysis' ...
            '\n[enter for default=', num2str(Th),']: '];
    hT = input(hTmess);
    if ~(isempty(hT))
         Th = hT;               % Th: Threshold for throwing away frames see line 34 in DDM_core2
    end    

     pentmess = ['Enter the start point for the characteristic time\n in the \I(tau) fits analysis [enter for default=',num2str(tempStPt),']:'];
     pent = input(pentmess);     % Show default
     if ~(isempty(pent))
          tempStPt = pent;               % Th: Threshold for throwing away frames see line 34 in DDM_core2
     end     

     querfmess = ['Enter the start point for the oscillation frequency\n in the \I(tau) fits analysis \n[enter for default=',num2str(freqStPt),']:'];
     querf = input(querfmess);     % Show default
     if ~(isempty(querf))
          freqStPt = querf;               % Th: Threshold for throwing away frames see line 34 in DDM_core2
     end    
     
     nufmess = ['Enter the number of the function to be fit for the \I(tau) fits analysis \n[enter for default=',num2str(tempStPt),']:'];
     nuf = input(nufmess);     % Show default
     if ~(isempty(nuf))
          fun = nuf;               % Th: Threshold for throwing away frames see line 34 in DDM_core2
     end    
     
    oofmess = ['Enter the start point for the \tau_q analysis \n[enter for default=',num2str(fixed),']:'];
    oof = input(oofmess);     % Show default
    if ~(isempty(oof))
          StPtfoo = oof;               % Th: Threshold for throwing away frames see line 34 in DDM_core2
    end    
     
    dexmess = ['Enter the fixed exponent for the \tau_q analysis \n[enter for default=',num2str(fixed),']:'];
    dexif = input(dexmess);     % Show default
    if ~(isempty(dexif))
          fixed = dexif;               % Th: Threshold for throwing away frames see line 34 in DDM_core2
    end    

end

% Check if inputs correct, if not revise interactively
while ((check~='y') && (check~='Y'))
    disp([' ']);
    disp(['- The subdirectory of MicroscopeData selected'])
        disp(['is "', userdir,'".']);
    disp(['- The microscope''s micron/px ratio = ', num2str(px2mum),'.']);
    disp(['- The video''s 1/framerate = ', num2str(frames2s),'.']);
    disp('- The threshold for throwing away')
        disp(['frames in the analysis = ', num2str(Th),'.']);
    disp('- The start point for the I(\tau) fits analysis')
        disp(['is ',num2str(tempStPt),'.']);
    disp('- The function to be fitted to the saturation tau');
        disp(['functions is function number ', num2str(fun)]);
    disp(['- The start point for the \tau_q fits is ',num2str(Stptfoo),'.'])
    disp('- The fixed exponent for the ')
        disp(['\tau_q analysis is ',num2str(fixed),'.']);
    
    while ((check ~= 'n') && (check ~= 'y') && (check ~= 'Y') && (check ~= 'N'))
             peck = input('\nIs this correct [y/n]? ', 's');
             if(isempty(peck))
                 check = 'y';
             else
                check = peck;
             end    
             disp([ check ]);
    end
       
    if (strcmp(check,'n')) || (strcmp(check, 'N'))
        while ((check ~= 'd') && (check ~= 'f') && (check ~= 'r') && (check ~= 't') ....
                && (check ~= 'D') && (check ~= 'F') && (check ~= 'R') && (check ~= 'T'))
            check = input('\nWhich is wrong [d, r, f, t, e]? ','s');
        end
        switch check
            case {'d', 'D'} 
                userdir = uigetdir('C:\MicroscopeData','Please select the folder that contains the desired microscope video,');    
            case {'r', 'R'}
                px2mum = input('Enter the microscope''s micron/px ratio: ');
            case {'f', 'F'}
                framerate = input('Enter the framerate of the video: ');
            case {'t', 'T'}
                Th = input('Enter the threshold for throwing away \nframes in the analysis: ');
            case {'e', 'E,'}
                fixed = input('Enter the fixed exponent for the \tau_q analysis: ');           
            otherwise
                disp('error')
        end
    end
end

if isempty(list)        % regardless of initial state of framerate,
   ext = 'mat';            % Save frames2s with video if have not already done so
   filename = strjoin({'vidpars',ext},'.');
   path = strjoin({userdir,filename}, filesep);
   variables = {'frames2s','px2mum'};
   save(path, variables{:});
end

%% Call External Functions
fs = ometiffreader2(userdir);                % read in image stack

[ DDMdir, Iqtau ] = DDM_core2(fs, userdir);  % do Fourier Transforms, save Iqtau to folder under userdir (DDMdir), with timing


fit_Iqtau_colloids2(fun, frames2s, Iqtau, DDMdir);    % Do fitting of I(tor(q)) \propto tor(q) over q 
                                                                                              % save to folder under DDMdir (fitdir)

end