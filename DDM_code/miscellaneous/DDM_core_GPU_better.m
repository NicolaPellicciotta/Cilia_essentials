function [ Iqtau, err_Iqtau ] = DDM_core_GPU_better( frame_stack )
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

frame_stack = single(frame_stack);


%% general-purpose variables

max_tau = floor(N_frames/2);    %also equal to the number of lags
max_q = floor(frame_size/2);

% distance map

jj = repmat(gpuArray(1:frame_size),frame_size,1);
ii=jj';
cc = max_q+1;
distance_map =fftshift(round(sqrt((ii-cc).*(ii-cc)+(jj-cc).*(jj-cc)))+1); %if we fftshit here is only 1 time instead of Nframes/2
distance_map = distance_map(:);

dist_counts = accumarray(distance_map,ones(frame_size*frame_size,1));
max_dist = length(dist_counts);
%% actual DDM calculation
%(difference of actual images, then |FFT|^2, average and azimuthal average)

Iqtau = zeros(max_q, max_tau, 'gpuArray');
err_Iqtau = zeros(max_q, max_tau, 'gpuArray');

for tau=1:max_tau
    
    dummy = gpuArray(zeros(max_dist, ceil((N_frames-tau)/tau) ,'single'));
    tc=0;
    for i=1:tau:N_frames-tau,
        tempdiff = fft2( frame_stack(:,:,i)-frame_stack(:,:,i+tau) ) * fft2_norm_factor;
        abs_FT_diff_image = real( tempdiff ).^2 + imag( tempdiff ).^2;
        tc = tc+1;
        dummy(:,tc) = accumarray(distance_map,abs_FT_diff_image(:));
    end
    
    dummy = bsxfun(@rdivide,dummy,dist_counts); %same as dummy = dummy .* repmat(dist_counts,1,tc);
    oneD_power_spectrum = mean(dummy,2);
    Iqtau(:,tau) = oneD_power_spectrum(2:max_q+1);
    
    if nargout > 1
        std_oneD_power_spectrum = std(dummy,1,2);
        err_Iqtau(:,tau) = std_oneD_power_spectrum(2:max_q+1);
    end
    
end

Iqtau = gather(Iqtau);
err_Iqtau = gather(Iqtau);

end



