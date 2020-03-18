% this version will generate a square lattice of oscillating colloids

close all
clear all
clearvars -except tumble_prob tpc tpv
tic
%% constants

TotalTime = 30; % in seconds
FrameRate = 100; %fps
N_frames = TotalTime * FrameRate;
dt = 1/FrameRate;

frame_size = 2048;
DDM_FOV = 512;
min_DDM_window = frame_size/2-DDM_FOV/2 +1;
max_DDM_window = frame_size/2+DDM_FOV/2;


colloid_diameter = 6;

N_colloids = 8192;

xc0 = ceil(frame_size*rand(N_colloids,1));
yc0 = ceil(frame_size*rand(N_colloids,1));

angle0 = 2*pi*rand(N_colloids,1);
vel = 50; %velocity, px/s

tumble_prob = 1;%4e-3; %probability of a tumble;

 %% generate colloid's template

[X,Y]=meshgrid(-16:16,-16:16);
colloid = double(uint8( 64* (1-tanh( (sqrt(X.^2+Y.^2) - colloid_diameter/2) / 2 ))./2 )) ;
% imagesc(colloid)
% shg

%% generate positions

time_vec = (1:N_frames)/FrameRate;

tumble_matrix = rand(N_colloids,N_frames) <= tumble_prob; % element i,j is 1 if the ith colloid tumbles at the jth frame
tumbling_angle = double(tumble_matrix);
tumbling_angle(tumble_matrix) = 2*pi*rand(sum(tumble_matrix(:)),1); % element i,j is the change in angle that the ith colloid makes between frame j-1 and j
angle = repmat(angle0,1,N_frames) + cumsum(tumbling_angle,2);
dXc = vel*cos(angle)*dt; % matrix of displacements
dYc = vel*sin(angle)*dt;

Xc = cumsum([xc0, dXc],2); %trajectory matrix
Yc = cumsum([yc0, dYc],2);

Xc = mod( round(Xc)-1, frame_size) +1;
Yc = mod( round(Yc)-1, frame_size) +1;

index_positions_in_DDM_window = ( Xc >= min_DDM_window & Xc <= max_DDM_window & Yc >= min_DDM_window & Yc <= max_DDM_window);

%%
Xc(~index_positions_in_DDM_window) = NaN;
Yc(~index_positions_in_DDM_window) = NaN;
%%
Xc = Xc - min_DDM_window +1;
Yc = Yc - min_DDM_window +1;

%% generate frame_stack by convolving positions and colloid template
toc
comp_time=uint32(1:N_frames);
frame_stack = arrayfun( @(i) uint8(conv2(full(sparse( Yc(~isnan(Yc(:,i)),i),Xc(~isnan(Yc(:,i)),i),ones(sum(~isnan(Yc(:,i))),1),DDM_FOV,DDM_FOV  )),colloid,'same')) ,comp_time,'UniformOutput',false);
% frame_stack = arrayfun( @(i) uint8(full(sparse( Yc(:,i),Xc(:,i),256*ones(N_colloids,1),frame_size,frame_size  ))) ,comp_time,'UniformOutput',false);
frame_stack = cat(3,frame_stack{:});
toc
%%

% frame_stack = frame_stack((frame_size-DDM_FOV)/2+1:(frame_size+DDM_FOV)/2,(frame_size-DDM_FOV)/2+1:(frame_size+DDM_FOV)/2,:);

%%

std_fs=std(single(frame_stack),1,3);
imagesc(std_fs); shg

%%
tic;[Iqtau,err_Iqtau]=DDM_core(frame_stack);toc
imagesc(Iqtau)

% save(['ballistic_colloids_v=',num2str(vel),'_s=',num2str(0),'_TumbleProb=',num2str(tumble_prob),'.mat']);
beep