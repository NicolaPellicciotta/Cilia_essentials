% this version will generate a square lattice of oscillating colloids

close all
clear all
tic
%% constants

TotalTime = 30; % in seconds
FrameRate = 100; %fps
N_frames = TotalTime * FrameRate;
dt = 1/FrameRate;

sample_size = 2048;
sample_depth = 512;
DDM_FOV = 512;
focus_depth = 31;
min_DDM_window = sample_size/2-DDM_FOV/2 +1;
max_DDM_window = sample_size/2+DDM_FOV/2;
min_focus_depth = round((sample_depth-focus_depth)/2+1);
max_focus_depth = round((sample_depth+focus_depth)/2);

xm = -7:7;
ym = -7:7;
zm = -15:15;
[Xm,Ym,Zm] = meshgrid(xm,ym,zm);
PSF = exp(-( (Xm/2).^2 + (Ym/2).^2 + (Zm/6).^2   ));

colloid_diameter = 6;

N_colloids = 65536;

xc0 = ceil(sample_size*rand(N_colloids,1));
yc0 = ceil(sample_size*rand(N_colloids,1));
zc0 = ceil(sample_depth*rand(N_colloids,1));


theta0 = 2*pi*rand(N_colloids,1);
phi0 = pi*rand(N_colloids,1);
vel = 50; %velocity, px/s

tumble_prob = 0.01;%4e-3; %probability of a tumble;

 %% generate colloid's template

[X,Y,Z]=meshgrid(-16:16,-16:16,-16:16);
colloid = false(size(X));
colloid(X.*X + Y.*Y + Z.*Z < (colloid_diameter/2)^2) = true;
% imagesc(colloid)
% shg

%% generate positions

time_vec = (1:N_frames)/FrameRate;

tumble_matrix = rand(N_colloids,N_frames) <= tumble_prob; % element i,j is 1 if the ith colloid tumbles at the jth frame

tumbling_theta = tumble_matrix;
tumbling_phi = tumble_matrix;

tumbling_theta(tumble_matrix) = 2*pi*rand(sum(tumble_matrix(:)),1); % element i,j is the change in angle that the ith colloid makes between frame j-1 and j
tumbling_phi(tumble_matrix) = 2*pi*rand(sum(tumble_matrix(:)),1);

theta = repmat(theta0,1,N_frames) + cumsum(tumbling_theta,2);
phi = repmat(phi0,1,N_frames) + cumsum(tumbling_phi,2);

dXc = vel*sin(phi).*cos(theta)*dt; % matrix of displacements
dYc = vel*sin(phi).*sin(theta)*dt;
dZc = vel*cos(phi)*dt;

Xc = cumsum([xc0, dXc],2);  clear dXc;%trajectory matrix
Yc = cumsum([yc0, dYc],2);  clear dYc;
Zc = cumsum([zc0, dZc],2);  clear dZc;

Xc = mod( round(Xc)-1, sample_size) +1;
Yc = mod( round(Yc)-1, sample_size) +1;
Zc = mod( round(Zc)-1, sample_depth) +1;

index_positions_in_sampling_window = ( Xc >= min_DDM_window & Xc <= max_DDM_window & Yc >= min_DDM_window & Yc <= max_DDM_window & Zc >= min_focus_depth & Zc <= max_focus_depth);

%%
Xc(~index_positions_in_sampling_window) = NaN;
Yc(~index_positions_in_sampling_window) = NaN;
Zc(~index_positions_in_sampling_window) = NaN;

%%
Xc = Xc - min_DDM_window +1;
Yc = Yc - min_DDM_window +1;
Zc = Zc - min_focus_depth +1;


%% generate frame_stack by convolving positions and colloid template
toc
tic
comp_time=(1:N_frames);
dummy_sample = false(DDM_FOV,DDM_FOV, focus_depth);
PSF=gpuArray(PSF);
% frame_stack = cell(length(comp_time));
frame_stack =...
    arrayfun(...
        @(i) uint16(...
        gather(...
            sum(...
                convn(...
                    gpuArray(single(imdilate(update_FOV_3d(dummy_sample,Xc(:,i),Yc(:,i),Zc(:,i)),colloid))),PSF,'same'...
                )...
            ,3)...
            )...
        ),...
    comp_time,'UniformOutput',false);
frame_stack = cat(3,frame_stack{:});


% frame_stack = arrayfun( @(i) uint8(full(sparse( Yc(:,i),Xc(:,i),256*ones(N_colloids,1),frame_size,frame_size  ))) ,comp_time,'UniformOutput',false);
toc


