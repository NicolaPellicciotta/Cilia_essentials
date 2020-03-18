function [ Iqtau, err_Iqtau ] = DDM_core_GPU( frame_stack )
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

%% actual DDM calculation
%(difference of actual images, then |FFT|^2, average and azimuthal average)

Iqtau = zeros(max_q, max_tau, 'gpuArray');
err_Iqtau = zeros(max_q, max_tau, 'gpuArray');

for tau=1:max_tau
    
    accum_abs_FT_diff_image = zeros(frame_size,'gpuArray');
       
    for i=1:tau:N_frames-tau,
        tempdiff = fft2( frame_stack(:,:,i)-frame_stack(:,:,i+tau) ) * fft2_norm_factor;
        accum_abs_FT_diff_image = accum_abs_FT_diff_image +...
            real( tempdiff ).^2 +...
            imag( tempdiff ).^2;
    end
      
    averaged_abs_FT_diff_image = accum_abs_FT_diff_image./ceil((N_frames-tau)/tau); 
    oneD_power_spectrum = accumarray(distance_map,averaged_abs_FT_diff_image(:))./dist_counts;
    
    if nargout == 2
        reconstructed_profile_matrix = oneD_power_spectrum(distance_map);
        square_difference_matrix = (averaged_abs_FT_diff_image(:)-reconstructed_profile_matrix).^2;
        dummy_std = sqrt(1./(dist_counts-1).*accumarray(distance_map,square_difference_matrix(:)));
        
        err_Iqtau(:,tau) = dummy_std(2:max_q+1);        
    end
    Iqtau(:,tau) = oneD_power_spectrum(2:max_q+1);
    
end

Iqtau = gather(Iqtau);
err_Iqtau = gather(Iqtau);

end



