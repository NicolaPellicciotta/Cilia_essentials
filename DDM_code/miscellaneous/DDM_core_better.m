function [ Iqtau, err_Iqtau ] = DDM_core_better( frame_stack )
%DDM_core Computational core for DDM calculation
%   Detailed explanation goes here

%% input check

if size(frame_stack,1)~=size(frame_stack,2)
    disp('Warning: non-square frames, only a square region will be actually analysed');
    frame_size = min([size(frame_stack,1),size(frame_stack,2)]);
    frame_stack = frame_stack(1:frame_size,1:frame_size,:);
end

frame_size = min([size(frame_stack,1),size(frame_stack,2)]);
N_px_frame = frame_size*frame_size;
fft2_norm_factor = 1/N_px_frame;
N_frames = size(frame_stack,3);

%% general-purpose variables

max_tau = floor(N_frames/2);    %also equal to the number of lags
max_q = floor(frame_size/2);

%% distance map for fast radial average

jj = repmat((1:frame_size),frame_size,1);
ii=jj';
cc = max_q+1;
distance_map =fftshift(round(sqrt((ii-cc).*(ii-cc)+(jj-cc).*(jj-cc)))+1); %if we fftshit here is only 1 time instead of Nframes/2
distance_map = distance_map(:);

dist_counts = accumarray(distance_map,ones(frame_size*frame_size,1));
max_dist = length(dist_counts);

%% actual DDM calculation
%(difference of frames, then |FFT|^2, average and radial average)

Iqtau = zeros(max_q, max_tau, 'double');
err_Iqtau = zeros(max_q, max_tau, 'double');

for tau=1:max_tau
    
    dummy = zeros(max_dist, ceil((N_frames-tau)/tau) ,'single');
    tc=0;
    for i=1:tau:N_frames-tau, % 1d fourier transform for each initial time, averaged just after the for loop
        tempdiff = fftn( single(frame_stack(:,:,i))-single(frame_stack(:,:,i+tau)) ) * fft2_norm_factor; %fft2 of difference image
        abs_FT_diff_image = real( tempdiff ).^2 + imag( tempdiff ).^2;
        tc = tc+1;
        dummy(:,tc) = accumarray(distance_map,abs_FT_diff_image(:));	%radial average
    end
    
    dummy = bsxfun(@rdivide,dummy,dist_counts); %same as dummy = dummy .* repmat(dist_counts,1,tc);
    oneD_power_spectrum = mean(dummy,2);
    Iqtau(:,tau) = oneD_power_spectrum(2:max_q+1);	%fill each column of the output with the 1D power spectrum just calculated
    
    if nargout > 1
        std_oneD_power_spectrum = std(dummy,1,2);
        err_Iqtau(:,tau) = std_oneD_power_spectrum(2:max_q+1);
    end
end
%{
%% actual DDM calculation
%(difference of frames, then |FFT|^2, average and radial average)

Iqtau = zeros(max_q, max_tau, 'double');
err_Iqtau = zeros(max_q, max_tau, 'double');

tempdiff = @(i,tau) = fftn( single(frame_stack(:,:,i))-single(frame_stack(:,:,i+tau)) ) * fft2_norm_factor;

for tau=1:max_tau
    
    dummy = zeros(max_dist, ceil((N_frames-tau)/tau) ,'single');
    tc=0;
    for i=1:tau:N_frames-tau, % on-the-fly average of fft2 of difference images
        tempdiff = fftn( single(frame_stack(:,:,i))-single(frame_stack(:,:,i+tau)) ) * fft2_norm_factor; %fft2 of difference image
        abs_FT_diff_image = real( tempdiff ).^2 + imag( tempdiff ).^2;
        tc = tc+1;
        dummy(:,tc) = accumarray(distance_map,abs_FT_diff_image(:));	%radial average
    end
    
    dummy = bsxfun(@times,dummy,dist_counts); %same as dummy = dummy .* repmat(dist_counts,1,tc);
    oneD_power_spectrum = mean(dummy,2);
    Iqtau(:,tau) = oneD_power_spectrum(2:max_q+1);	%fill each column of the output with the 1D power spectrum just calculated
    
    if nargout > 1
        std_oneD_power_spectrum = std(dummy,1,2);
        err_Iqtau(:,tau) = std_oneD_power_spectrum(2:max_q+1);
    end
end
%}
end



