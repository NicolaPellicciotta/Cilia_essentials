% this version will generate a random lattice of oscillating dots

close all
clear
tic

ampl_vec = [1:1:60];
% ampl_vec = [1,60];

for sc = numel(ampl_vec):-1:1
    
%% constants

Sim(sc).TotalTime = 10; % in seconds
Sim(sc).FrameRate = 20; %fps
Sim(sc).N_frames = Sim(sc).TotalTime * Sim(sc).FrameRate;

Sim(sc).frame_size = 2048;
Sim(sc).border = 32;

dummy_frame_size = Sim(sc).frame_size+2*Sim(sc).border;

% nndist = 32; %distance between next neighbours
% xc0 = nndist/2 : nndist : dummy_frame_size;
% yc0 = nndist/2 : nndist : dummy_frame_size;
% N_oscillators = length(xc0)^2;


Sim(sc).N_oscillators = 8192;

Sim(sc).xc0 = ceil(dummy_frame_size*rand(Sim(sc).N_oscillators,1)); %random positions
Sim(sc).yc0 = ceil(dummy_frame_size*rand(Sim(sc).N_oscillators,1));
[Sim(sc).xc0, sort_index] = sort(Sim(sc).xc0);
Sim(sc).yc0 = Sim(sc).yc0(sort_index);


Sim(sc).oscillator_diameter = 12;
% Sim(sc).osc_ampl = Sim(sc).oscillator_diameter .* ones(Sim(sc).N_oscillators,1) .* 2;
Sim(sc).osc_ampl = ampl_vec(sc) .* ones(Sim(sc).N_oscillators,1);

Sim(sc).f = 5 ;%+ 0.2*randn(Sim(sc).N_oscillators,1);  %mean frequency (Hz)
% Sim(sc).angle = 2*pi*rand(Sim(sc).N_oscillators,1);
% Sim(sc).angle = pi/6 + pi/32*randn(sqrt(Sim(sc).N_oscillators));
Sim(sc).angle = 0; %in radians
% Sim(sc).phase = 8*pi/Sim(sc).frame_size .* Sim(sc).xc0; %phase has to be defined wrt the real size of the video % this is metachronal wave
Sim(sc).phase = 2*pi*rand(Sim(sc).N_oscillators,1);   %no metachronal wave

%% generate oscillator's template

[X,Y]=meshgrid(-16:16,-16:16);
Sim(sc).template = double(uint8( 20* (1-tanh( (sqrt(X.^2+Y.^2) - Sim(sc).oscillator_diameter/2) / 3 ))./2 )) ;
% imagesc(colloid)
% shg

%% generate positions

Sim(sc).time_vec = (1:Sim(sc).N_frames)/Sim(sc).FrameRate;

% Sim(sc).X_pos_fun = @(time_vec) mod(round(Sim(sc).xc0 + cos(Sim(sc).angle).*Sim(sc).osc_ampl.*     ( sin(2*pi*Sim(sc).f.*time_vec + Sim(sc).phase) + 0.1.*sin(2* (2*pi*Sim(sc).f.*time_vec + Sim(sc).phase)) )     )-1,dummy_frame_size)+1;%asymmetric waveform
% Sim(sc).Y_pos_fun = @(time_vec) mod(round(Sim(sc).yc0 + sin(Sim(sc).angle).*Sim(sc).osc_ampl.*     ( sin(2*pi*Sim(sc).f.*time_vec + Sim(sc).phase) + 0.1.*sin(2* (2*pi*Sim(sc).f.*time_vec + Sim(sc).phase)) )     )-1,dummy_frame_size)+1;

Sim(sc).X_pos_fun = @(time_vec) mod(round(Sim(sc).xc0 + cos(Sim(sc).angle).*Sim(sc).osc_ampl.*      sin(2*pi*Sim(sc).f.*time_vec + Sim(sc).phase)      )-1,dummy_frame_size)+1;
Sim(sc).Y_pos_fun = @(time_vec) mod(round(Sim(sc).yc0 + sin(Sim(sc).angle).*Sim(sc).osc_ampl.*      sin(2*pi*Sim(sc).f.*time_vec + Sim(sc).phase)      )-1,dummy_frame_size)+1;

Xc = arrayfun(Sim(sc).X_pos_fun,  Sim(sc).time_vec,'UniformOutput',false);
Yc = arrayfun(Sim(sc).Y_pos_fun,  Sim(sc).time_vec,'UniformOutput',false);

%%
toc
frame_stack = arrayfun( @(Yc,Xc) uint8(conv2(full(sparse( Yc{:},Xc{:},ones(Sim(sc).N_oscillators,1),dummy_frame_size,dummy_frame_size  )),Sim(sc).template,'same')) ,Yc, Xc,'UniformOutput',false);
frame_stack = cat(3,frame_stack{:});
toc

frame_stack = frame_stack(Sim(sc).border+1:Sim(sc).frame_size+Sim(sc).border,Sim(sc).border+1:Sim(sc).frame_size+Sim(sc).border,:);
toc

Sim(sc).Iqtau = DDM_core(frame_stack);
figure(sc)
imagesc(Sim(sc).Iqtau); colorbar;title(['Amplitude = ',num2str(ampl_vec(sc)),' px']);
% Sim(sc).Iqvectau = DDM_core_anisotropy(frame_stack,false,[],30);
toc

end
save('changing_amplitude_oscillation_2','Sim','-v7.3');

%% for investigating amplitude
close all
load('changing_amplitude_oscillation_2.mat')
%%
figure(10);
set(gcf,'Units','Normalized','Position',[0 0 1 1])

for sc=1:numel(Sim)
    clf
    Sim(sc).Iq = mean(Sim(sc).Iqtau,2);
    plot(Sim(sc).Iq)
    hold on
    smoothed_Iq = smooth(Sim(sc).Iq,15);
    findpeaks(smoothed_Iq);
    [pk,du] = findpeaks(smoothed_Iq);
    mp(sc) = du(pk == max(pk));
    drawnow
    pause(0.1);
end

Iq_vs_ampl = horzcat(Sim(:).Iq);

if ~exist('ampl_vec','var')
    ampl_vec = mean(horzcat(Sim(:).osc_ampl));
end

mode_vec = [1:512]';
qVec = 2*pi/1024.*mode_vec;
lambda_vec = 1024./mode_vec;

figure(11)
plot(ampl_vec,mp,'o')
ylabel('mode of DDM peak')


figure(12)
plot(ampl_vec,1024./mp,'o')
ylabel('lambda of DDM peak')

cftool(2*ampl_vec,1024./mp)


% figure(1)
% plot3(repmat(qVec,1,numel(ampl_vec)),repmat(ampl_vec,numel(qVec),1),Iq_vs_ampl)
% xlim([0,1])
% 
% figure(2)
% plot3(repmat(mode_vec,1,numel(ampl_vec)),repmat(ampl_vec,numel(qVec),1),Iq_vs_ampl)
% xlim([0, 100])
% 
% figure(3)
% plot3(repmat(lambda_vec,1,numel(ampl_vec)),repmat(ampl_vec,numel(qVec),1),Iq_vs_ampl)
% xlim([0 400])