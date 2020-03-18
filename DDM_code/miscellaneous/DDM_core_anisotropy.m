function [ Iqtau, err_Iqtau ] = DDM_core_anisotropy( frame_stack, flag_GPU, max_tau, angle_bins_width)
 
%DDM_core_anisotropy Computational core for DDM calculation
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


if nargin < 4
    angle_bins_width = 1;
end
if nargin < 3 || isempty(max_tau)
    max_tau = floor(N_frames/2);
end
if nargin < 2
    flag_GPU = false;
end



%% general-purpose variables

max_q = floor(frame_size/2);

% distance map

jj = repmat([1:frame_size],frame_size,1);
ii=jj';
cc = max_q+1;

distance_map = fftshift(round(sqrt((ii-cc).*(ii-cc)+(jj-cc).*(jj-cc)))+1); %if we fftshit here is only 1 time instead of Nframes/2
distance_map = distance_map(:);

angle_map = fftshift( round( 1/angle_bins_width * atand( (ii-cc)./(jj-cc) )) * angle_bins_width );
angle_map = mod(angle_map(:), 180)+1;
angle_map(1,1) = min(angle_map(:));

counts = accumarray({distance_map, angle_map(:)},ones(frame_size*frame_size,1));

unique_angles = unique(angle_map);
n_angle_bins = length(unique_angles);

%% actual DDM calculation
%(difference of actual images, then |FFT|^2, average and azimuthal average)

Iqtau = zeros(max_q, n_angle_bins, max_tau, 'double');
err_Iqtau = zeros(max_q, n_angle_bins, max_tau, 'double');

for tau=1:max_tau
    if flag_GPU
        accum_abs_FT_diff_image = zeros(frame_size,'gpuArray');
    else
        accum_abs_FT_diff_image = zeros(frame_size,'double');
    end
    for i=1:tau:N_frames-tau, 
        tempdiff = fft2( double(frame_stack(:,:,i))-double(frame_stack(:,:,i+tau)) ) * fft2_norm_factor;
        accum_abs_FT_diff_image = accum_abs_FT_diff_image +...
            real( tempdiff ).^2 +...
            imag( tempdiff ).^2;
    end

    
    
    if flag_GPU
        dummy = gather(accum_abs_FT_diff_image./ ceil((N_frames-tau)/tau));
    else
        dummy = accum_abs_FT_diff_image./ ceil((N_frames-tau)/tau);
    end
    dummy_polar = accumarray({distance_map, angle_map},dummy(:))./counts;
    
    
    if nargout == 2
        reconstructed_profile_matrix = dummy_polar(sub2ind(size(dummy_polar),distance_map,angle_map));
        square_difference_matrix = (dummy(:)-reconstructed_profile_matrix).^2;
        dummy_std = sqrt(1./(counts-1).*accumarray({distance_map, angle_map},square_difference_matrix(:)));
        
        err_Iqtau(:,:,tau) = dummy_std(2:max_q+1,unique_angles);
        
    end
    Iqtau(:,:,tau) = dummy_polar(2:max_q+1,unique_angles);
    

    % the C/C++ style average is needed to avoid building a big
    % temporary matrix we would have to average on in any case
end

end





