


mo=moviereader('40X_X100Y275.04Aug2017_16.58.30.movie');
%%
clear all;
bsz = 4;
load('changing_freq_distr_q_01_fs.mat')
for sc=1:numel(Sim)
    
fs=Sim(sc).fs;
fs=fs(:,:,1:1000);
% parameters
height = size(fs,1);
width = size(fs,2);
Nframes = size(fs,3);
FR = Sim(1).FrameRate%  mo.FrameRate; % frame rate

% fft window
window = hann(floor(Nframes/2));
% wwindow = hann(floor(Nframes/2));

% find how long will the frequency vector be
[~,dummyf] = periodogram(AutoCorr(double(squeeze(fs(1,1,:))),floor(Nframes/2)),window,floor(Nframes/2),FR);
Nfreqs = numel(dummyf);
clear dummyf;

% create binning map

binningmap = create_binning_map([height, width],bsz);

% N_frames, height*box matrix
freq_ind = repmat((1:Nfreqs)',1,height*bsz);

% row vector with the index within a binning coloumn (always start from 1,
% then proper placement in frequency map will be done by frequency_map(:,cc) = frequency_strip;
temp_ind = reshape(binningmap(:,1:bsz),height * bsz,1)';
% put into matrix form for accumarray
temp_ind_mat = repmat(temp_ind,Nfreqs,1);
    
    
% initialise frequency map
frequency_map = nan(height/bsz,width/bsz);

% for loop on columns of blocks

for cc = 1:ceil(width/bsz)
    fprintf('%.2d/%.2d', cc, ceil(width/bsz) );
    ccleft = (cc-1) * bsz + 1;
    ccright = ccleft + bsz-1;
    
    % take the first block, reshape it to be a N_frames, height*box matrix
    % (good for periodogram)
    % Each column is a time signal
    temp_fs = reshape(fs(:,ccleft:ccright,:),height * bsz,Nframes)'; 
    
    % subtract the mean
    temp_fs = double(temp_fs) - mean(temp_fs,'double');
    
    % take autocorrelation over time
    temp_fs = AutoCorr(temp_fs,floor(Nframes/2));
    
    % make Fourier transforms
    % pxx is a matrix with the spectrum of each pixel as a column
    [pxx, frequencies] = periodogram(temp_fs,window,size(temp_fs,1),FR);
%     [pwxx, frequencies] = pwelch(temp_fs,wwindow,[],numel(wwindow),FR);

    % average points of the spectra with same box index and same frequency
    boxavg_pxx = accumarray({freq_ind(:),temp_ind_mat(:)},pxx(:)) ./ (bsz*bsz);
    
    % times a downsampled spectrum
    % resample by averaging two close frequencies
%     resampled_boxavg_pxx = bin_matrix(boxavg_pxx,[2 1]); % maybe easier to do boxavg_pxx(1:2:end,:)
%     resampled_frequencies = bin_matrix(frequencies,[2 1]);
%     prod_boxavg_pxx = boxavg_pxx(1:size(resampled_boxavg_pxx,1),:) .* resampled_boxavg_pxx;
%     prod_frequencies = frequencies(1:size(resampled_boxavg_pxx,1));
    
    % peak detection
    frequency_strip = nan(size(boxavg_pxx,2),1);
    parfor ii = 1:size(boxavg_pxx,2)
        if any(isnan(abs(boxavg_pxx(:,ii))))==0;
        % find peaks
        [pks,locs] = findpeaks(boxavg_pxx(:,ii),frequencies);
        
        % index higher peak
        [~,where] = max(pks);
        
        % write
        frequency_strip(ii) = locs(where);
        end
    end %parfor
    
    frequency_map(:,cc) = frequency_strip;
    fprintf(repmat('\b',1,5))
end
Sim(sc).frequency_map=frequency_map;

end

save('changing_freq_distr_q_01_fs_fft.mat');

%%