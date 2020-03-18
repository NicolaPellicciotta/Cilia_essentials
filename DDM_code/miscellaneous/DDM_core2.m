function [  DDMdir, Iqtau ] = DDM_core2( frame_stack, userdir )
%DDM_core Computational core for DDM calculation
%   Detailed explanation goes here

%%
disp('_ _ _ _ ');
disp('DDM_core2');
tic;                            % Start timing

%% Constants
frame_size = min([size(frame_stack,1),size(frame_stack,2)]);

params = csvread('Settings.txt',1,0); % px2mum, Th, tempStpt, freqStpt, fun, Stptfoo, fixed
Th = params(2);

%input check
if size(frame_stack,1)~=size(frame_stack,2)
    disp('Warning: non-square frames.')
    disp('Only a square region will be analysed');
    frame_stack = frame_stack(1:frame_size,1:frame_size,:);
end

N_px_frame = frame_size*frame_size;
fft2_norm_factor = 1/N_px_frame;
N_frames = size(frame_stack,3);
max_tau = floor(N_frames/2);    %also equal to the number of lags
max_q = floor(frame_size/2);

%% distance map for fast radial average
jj = repmat((1:frame_size),frame_size,1);
ii=jj';
cc = max_q+1;
distance_map = fftshift(round(sqrt((ii-cc).*(ii-cc)+(jj-cc).*(jj-cc)))+1); %if we fftshit here is only 1 time instead of Nframes/2
distance_map = distance_map(:);

dist_counts = accumarray(distance_map,ones(frame_size*frame_size,1));

%% actual DDM calculation
%(difference of frames, then |FFT|^2, average and radial average)
Iqtau = zeros(max_q, max_tau, 'double');   % Th = 100;   % Threshold for throwing away frames

for tau=1:max_tau
    ind_frames_1 = randi(tau,1):tau:N_frames-tau;
    ind_frames_2 = randi(N_frames-tau,[1,Th]);
    if numel(ind_frames_2)<numel(ind_frames_1)
        ind_frames = ind_frames_2;
    else
        ind_frames = ind_frames_1;
    end
    
    temp_diff = single(frame_stack(:,:,ind_frames)) - single(frame_stack(:,:,ind_frames+tau));
    FT_temp_diff = fft2(temp_diff) * fft2_norm_factor;
    averaged_abs_FT_diff_image = mean(real(FT_temp_diff).^2 + imag(FT_temp_diff).^2, 3);
    oneD_power_spectrum = accumarray(distance_map,averaged_abs_FT_diff_image(:))./dist_counts;	%radial average
    Iqtau(:,tau) = oneD_power_spectrum(2:max_q+1);	%fill each column of the output with the 1D power spectrum just calculated
    
end
toc          % Stop timing, print

%% Save Iqtau matrix
% Save under
% userdir\DDM_core2_[Th]_[mmdd_HHMM]\DDM_core2_[Th]_[mmdd_HHMM].mat

t = datestr(now, 'mmdd_HHMM');
DDMname = ['DDM_core2_',num2str(Th),'_',t];
DDMdir = strjoin({userdir,DDMname}, filesep);
if ~isdir(DDMdir)
    mkdir(DDMdir);
end

ext = 'mat';
filename = strjoin({DDMname,ext},'.');
path = strjoin({DDMdir,filename}, filesep);
variables = {'Iqtau'};

save(path, variables{:});

end