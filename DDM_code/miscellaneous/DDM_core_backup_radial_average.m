function [ Iqtau, err_Iqtau ] = DDM_core( frame_stack, flag_GPU )
%DDM_core Computational core for DDM calculation
%   Detailed explanation goes here

%% input check

if size(frame_stack,1)~=size(frame_stack,2)
    disp('Warning: non-square frames, only a square region will be actually analysed');
    frame_size = min([size(frame_stack,1),size(frame_stack,2)]);
    frame_stack = frame_stack(1:frame_size,1:frame_size,:);
end
if nargin < 2
    flag_GPU = false;
end

frame_size = min([size(frame_stack,1),size(frame_stack,2)]);
N_px_frame = frame_size*frame_size;
fft2_norm_factor = 1/N_px_frame;
N_frames = size(frame_stack,3);

%% general-purpose variables

max_tau = floor(N_frames/2);    %also equal to the number of lags
max_q = floor(frame_size/2);

% distance map

jj = repmat([1:frame_size],frame_size,1);
ii=jj';
cc = max_q+1;
distance_map =fftshift(round(sqrt((ii-cc).*(ii-cc)+(jj-cc).*(jj-cc)))+1); %if we fftshit here is only 1 time instead of Nframes/2
distance_map = distance_map(:);

dist_counts = accumarray(distance_map,ones(frame_size*frame_size,1));

%% actual DDM calculation
%(difference of actual images, then |FFT|^2, average and azimuthal average)

Iqtau = zeros(max_q, max_tau, 'double');
err_Iqtau = zeros(max_q, max_tau, 'double');

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
    dummy_profile = accumarray(distance_map,dummy(:))./dist_counts;
    %store(tau)=mean(dummy(:));
    %     dummy_std = accumarray(distance_map,dummy(:),[],@std);
    
    if nargout == 2
        reconstructed_profile_matrix = dummy_profile(distance_map);
        square_difference_matrix = (dummy(:)-reconstructed_profile_matrix).^2;
        dummy_std = sqrt(1./(dist_counts-1).*accumarray(distance_map,square_difference_matrix(:)));
        
        err_Iqtau(:,tau) = dummy_std(2:max_q+1);
        
    end
    Iqtau(:,tau) = dummy_profile(2:max_q+1);
    
    %         Iqtau(:,tau) = azaver(gather(accum_abs_FT_diff_image./ ceil((N_frames-tau)/tau)), distance_map);
    %         Iqtau(:,tau) = azaver( accum_abs_FT_diff_image./ ceil((N_frames-tau)/tau) , distance_map);
    % the C/C++ style average is needed to avoid building a big
    % temporary matrix we would have to average on in any case
end

end

%% other versions of the code (generally slower or more RAM hungry)

%{

for tau=1:max_tau,
        Iqtau(:,tau) = azaver(mean(abs(fftshift(my_fft2(diff(double(frame_stack(:,:,1:tau:end)),1,3)))).^2,3));
        % a tad heavy on RAM beacause it makes a temporary matrix
        % double(diff(frame_stack(:,:,1:tau:end))), which is as large as
        % the original frame_stack when tau==1. However it is 1.5x faster
        % than the following RAM friendly version.
    end

%%%% 1
% Fourier transform of all images, very heavy on memory
FT_frame_stack = fftshift(fft2(frame_stack));  %fft2 on a 3d matrix is faster than the usual fft2, as it calls size(matrix,3) times fft2 on each slice, while the
Iqtau = zeros(max_q, max_tau, 'double');
for tau=1:max_tau         %this bit is similar to the multitau algorithm in MSD calculation
    Iqtau(:,tau) = mean( azaver( abs( diff(FT_frame_stack(:,:,1:tau:end),1,3) ).^2 ) );
end

%%%% 1bis
%(less RAM, but nested for-loop)
FT_frame_stack = fftshift(fft2(frame_stack));
Iqtau = zeros(max_q, max_tau, 'double');
for tau=1:max_tau         %this bit is similar to the multitau algorithm in MSD calculation
    for i=1:tau:N_frames-tau
        Iqtau(:,tau) = Iqtau(:,tau) + azaver( abs(FT_frame_stack(:,:,i) - FT_frame_stack(:,:,i+tau)).^2 );    % I could have used a temporary matrix, but this way we save a bit on memory
    end
    Iqtau(:,tau) = Iqtau(:,tau)./ceil((N_frames-tau)/tau);
end

%%%% 2
% this time I average over the differential images, and then azimuthally.
% it's the fastest, but also exceptionally heavy on RAM
FT_frame_stack = fftshift(fft2(frame_stack));
Iqtau = zeros(max_q, max_tau, 'double');
for tau=1:max_tau         %this bit is similar to the multitau algorithm in MSD calculation
    Iqtau(:,tau) = azaver( mean(abs(diff(FT_frame_stack(:,:,1:tau:end),1,3)).^2,3) );
end

%%%% 2bis
% variation of 2 (nested for-loop)
FT_frame_stack = fftshift(fft2(frame_stack));
Iqtau = zeros(max_q, max_tau, 'double');
for tau=1:max_tau         %this bit is similar to the multitau algorithm in MSD calculation
    abs_diff_FT_image = zeros(frame_size,'double');
    for i=1:tau:N_frames-tau
        abs_diff_FT_image = abs_diff_FT_image + abs(FT_frame_stack(:,:,i) - FT_frame_stack(:,:,i+tau)).^2;
    end
    Iqtau(:,tau) = azaver(abs_diff_FT_image/ceil((N_frames-tau)/tau));
end

%%%% 4
% difference of actual images, then |FFT|^2, azimuthal average and average on the differential images
% really slow and RAM hungry
Iqtau = zeros(max_q, max_tau, 'double');
for tau=1:max_tau
    Iqtau(:,tau) = mean(azaver(abs(fftshift(fft2(diff(double(frame_stack(:,:,1:tau:end)),1,3)))).^2));
end

%%%% 4bis
% variation of the 4 (nested for-loops)
Iqtau = zeros(max_q, max_tau, 'double');
for tau = 1:max_tau
    for i=1:tau:N_frames-tau
        Iqtau(:,tau) = Iqtau(:,tau) + azaver( abs(fftshift(fft2( double(frame_stack(:,:,i))-double(frame_stack(:,:,i+tau)) ))).^2 );
    end
    Iqtau(:,tau) = Iqtau(:,tau)./ceil((N_frames-tau)/tau);
end

%}




