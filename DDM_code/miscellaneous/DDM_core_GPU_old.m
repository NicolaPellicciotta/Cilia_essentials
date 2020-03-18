function [ Iqtau ] = DDM_core_GPU( movie_object )
%DDM_core Computational core for DDM calculation
%   Detailed explanation goes here

%% input check

dummy = movie_object.read(1);
first_frame = dummy.IM;

if size(first_frame,1)~=size(first_frame,2)
    disp('Warning: non-square frames, only a square region will be actually analysed');
    flag_rectangular = true;
    frame_size = min(size(first_frame));
else
    flag_rectangular = false;
    frame_size = size(first_frame,1);
end


if frame_size ~= 2048,
frame_size = 1024;%min([size(first_frame,1),size(first_frame,2)]);
end
N_frames = movie_object.NumberOfFrames;

%% general-purpose variables

max_tau = floor(N_frames/2);    %also equal to the number of lags
max_q = floor(frame_size/2);


%% actual DDM calculation
%(difference of actual images, then |FFT|^2, average and azimuthal average)

Iqtau = zeros(max_q, max_tau, 'double');
ith_frame = zeros(frame_size,'gpuArray');
iplustauth_frame = zeros(frame_size,'gpuArray');

for tau=1:max_tau
    accum_abs_FT_diff_image = zeros(frame_size,'gpuArray');
    for i=1:tau:N_frames-tau,
        dummy = movie_object.read(i);
        if flag_rectangular, dummyim = dummy.IM; ith_frame(:) = dummyim(1:frame_size,1:frame_size);
        else ith_frame(:) = dummy.IM;
        end
        
        dummy = movie_object.read(i+tau);
        if flag_rectangular, dummyim = dummy.IM; iplustauth_frame(:) = dummyim(1:frame_size,1:frame_size);
        else iplustauth_frame(:) = dummy.IM;
        end
        accum_abs_FT_diff_image = accum_abs_FT_diff_image + abs( fftshift(my_fft2( double(ith_frame)-double(iplustauth_frame) )) ).^2;
    end
    
    Iqtau(:,tau) = azaver(gather(accum_abs_FT_diff_image./ ceil((N_frames-tau)/tau)));
    
    % the C/C++ style average is needed to avoid building a big
    % temporary matrix we would have to average on in any case
end



end

%% other versions of the code (generally slower or more RAM hungry)

%{

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




