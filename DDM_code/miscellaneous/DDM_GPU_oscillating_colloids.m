% this version will compute the DDM of a generated on-the-run square lattice of oscillating colloids

close all
clear all

%% constants

TotalTime = 200; % in seconds
FrameRate = 300; %fps
N_frames = TotalTime * FrameRate;
max_tau = N_frames/2; %max lag time, in frames

frame_size = 128;
nndist = 32; %distance between next neighbours
colloid_diameter = 6;

osc_ampl = colloid_diameter * 2;

f = 15;  %mean frequency (Hz)

%% generate colloid's template

[X,Y]=meshgrid(-10:10,-10:10);
colloid = gpuArray( (1-tanh( (sqrt(X.^2+Y.^2) - colloid_diameter/2) / 1 ))./2 ) ;

%% distance map

% jj = repmat([1:frame_size],frame_size,1);
% ii=jj';
% cc = frame_size/2+1;
% distance_map =round(sqrt((ii-cc).*(ii-cc)+(jj-cc).*(jj-cc)))+1; %if we fftshit here is only 1 time instead of Nframes/2
% distance_map = distance_map(:);
% dist_counts = accumarray(distance_map,ones(frame_size*frame_size,1));

%% DDM

Iqtau = zeros(frame_size,max_tau,'gpuArray');

for tau = 1:max_tau %lag time, in frames
    
    % prepare the accumulative image
    accum_abs_FT_diff_image = zeros(frame_size,'gpuArray');
    
    % prepare the first frame
    positions_a = zeros(frame_size,'gpuArray');
    xa = gpuArray( osc_ampl/2 + (osc_ampl * sin(2*pi*f* 1/FrameRate) : nndist:frame_size));
    ya = gpuArray( osc_ampl/2 + (osc_ampl * sin(2*pi*f* 1/FrameRate) : nndist:frame_size));
    [Xa,Ya] = meshgrid(xa,ya);
    good_a = (Ya > 0 & Ya <= frame_size & Xa > 0 & Xa<= frame_size);
    positions_a(  sub2ind( [frame_size,frame_size], round(Ya(good_a)),round(Xa(good_a)) )  )=1;
    frame_a = conv2(positions_a,colloid, 'same');
    
    for tc = 1:tau:N_frames-tau %loop on the frames
        
        tb = (tc+tau)/FrameRate;
        %prepare the frame at tc+tau
        positions_b = zeros(frame_size,'gpuArray');
        xb = gpuArray( osc_ampl/2 +  (osc_ampl * sin(2*pi*f* tb) : nndist:frame_size) );
        yb = gpuArray( osc_ampl/2 +  (osc_ampl * sin(2*pi*f* tb) : nndist:frame_size) );
        [Xb,Yb] = meshgrid(xb,yb);
        good_b = (Yb > 0 & Yb <= frame_size & Xb > 0 & Xb<= frame_size);
    positions_b(  sub2ind( [frame_size,frame_size], round(Yb(good_b)),round(Xb(good_b)) )  )=1;
        frame_b = conv2(positions_b,colloid, 'same');
        
        tempdiff = 0.01*fft2( frame_a - frame_b );                       % fft of the differential image
        accum_abs_FT_diff_image = accum_abs_FT_diff_image +...      %let's accumulate it
            real( tempdiff ).^2 +...
            imag( tempdiff ).^2;
        
        frame_a = frame_b; %so I don't have to recompute frame_a
        
    end
        
    accum_abs_FT_diff_image = accum_abs_FT_diff_image./ ceil((N_frames-tau)/tau); %divide by the number of fft that I accumulated
    
    horizontal_profile = mean(accum_abs_FT_diff_image);
    Iqtau(:,tau) = horizontal_profile;
    
    
end

