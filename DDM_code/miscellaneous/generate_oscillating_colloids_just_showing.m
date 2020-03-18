% this version will generate a square lattice of oscillating colloids

close all
clear all

%% constants

TotalTime = 10; % in seconds
FrameRate = 300; %fps
N_frames = TotalTime * FrameRate;

frame_size = 1024;
nndist = 32; %distance between next neighbours
colloid_diameter = 6;

osc_ampl = colloid_diameter * 2;

f = 15;  %mean frequency (Hz)

%% generate colloid's template

[X,Y]=meshgrid(-10:10,-10:10);
colloid = uint8( 255* (1-tanh( (sqrt(X.^2+Y.^2) - colloid_diameter/2) / 1 ))./2 ) ;
% imagesc(colloid)
% shg

%% generate lattice

positions = false(frame_size);
[Xc,Yc] = meshgrid(osc_ampl/2:nndist:frame_size,osc_ampl/2:nndist:frame_size);
positions(Yc,Xc)=true;

%% generate first_frame

first_frame = uint8(conv2(single(positions),single(colloid), 'same'));
figure(1)
imshow(first_frame,[])

%% distance map

jj = repmat([1:frame_size],frame_size,1);
ii=jj';
cc = frame_size/2+1;
distance_map =round(sqrt((ii-cc).*(ii-cc)+(jj-cc).*(jj-cc)))+1; %if we fftshit here is only 1 time instead of Nframes/2
distance_map = distance_map(:);
dist_counts = accumarray(distance_map,ones(frame_size*frame_size,1));

%% update lattice, differential image and plots

centred_resized_figure(2,[2,2]);

for tc = 1:N_frames
    
    t = tc/FrameRate; %real time (sec)
    
    % update lattice
    current_positions = false(frame_size);
    [Xc,Yc] = meshgrid(osc_ampl/2 + round( osc_ampl * sin(2*pi*f*t)):nndist:frame_size,...
        osc_ampl/2 + round( osc_ampl * sin(2*pi*f*t)):nndist:frame_size);
    current_positions(Yc(Yc>0),Xc(Xc>0))=true;
    
    % update frame
    current_frame = uint8(conv2(single(current_positions),single(colloid), 'same'));
    
    % differential frame
    frame_diff = double(first_frame) - double(current_frame);
    
    % fourier transform
    fft_frame_diff = fftshift(fft2( frame_diff ));
    scattering_image = real(fft_frame_diff).^2 + imag(fft_frame_diff).^2;
    
    % plot
    subplot(221)
    subimage(mat2gray(frame_diff, [-255 255]));
    xlabel(['first frame - current frame, time = ',sprintf('%.3f',(t)),'s']);
    
    subplot(222)
    colormap('default')
    subimage(255 * mat2gray(scattering_image),colormap);
    xlabel('Fourier Transform of the differential image');
    
    subplot(223)
    plot(mean(scattering_image));
    xlabel('Average on the columns of the FT');
    
    % radial average plot
    subplot(224)
    plot(accumarray(distance_map,scattering_image(:)) ./ dist_counts);
    xlabel('Radial Average of the FT');
    
    pause(0.05)
end

