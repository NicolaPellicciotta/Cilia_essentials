clear all
close all


%% Constants
Size = 480;
SizeCell = 80;
FOVsizevect = [4 8 10 16 20 24 32 40 48 60 80 96 120 160 240 480];

% Size = 1024;
% SizeCell = 26; % 2 times the side if hexagonal patch
% FOVsizevect = [8 16 32 64 128 256 512 1024];

TotalTime = 10; % in seconds
FrameRate = 300; %fps
N_frames = TotalTime*FrameRate; %total number of frames

f = 15;  %mean frequency
sigma = 0.2; %sigma

%% Create hexagonal patches

CellRadius = SizeCell/2;

Y=round(1:(CellRadius/2*sqrt(3)):Size+SizeCell);
X=repmat([1:3*CellRadius:Size+SizeCell;3/2*CellRadius+1:3*CellRadius:Size+SizeCell],ceil(length(Y)/2),1);
if size(X,1) > length(Y), X(end,:)=[]; end
Y=repmat(Y,size(X,2),1);
Y=Y';
A=zeros(Size+SizeCell); A(sub2ind([Size+SizeCell,Size+SizeCell],round(Y(:)),round(X(:))))=1;imagesc(A);axis image
centers=[X(:),Y(:)];
[x,y]=meshgrid((1:Size+SizeCell),(1:Size+SizeCell));
xy=[x(:),y(:)];
blbl=knnsearch(centers,xy);
indices = reshape(blbl,Size+SizeCell,Size+SizeCell);
indices = indices(1:Size,1:Size);

%% Parameters initialisation (for hexagonal patches)

AmplMat = 127*ones(Size); %for uint8 images
CIMat = 128*ones(Size);

PhaseMat = 2*pi*rand(Size); %initial phase in angle units

N_Cells = max(indices(:));
Freq = f+sigma*randn(N_Cells,1);
FreqMat = Freq(indices);
FreqMat=reshape(FreqMat,Size,Size);
%}
%% Parameters initialisation (for square patches)
%{
AmplMat = 127*ones(Size); %for uint8 images
CIMat = 128*ones(Size);

PhaseMat = 2*pi*rand(Size); %initial phase in angle units


Freq = f+sigma*randn(Size/SizeCell); %matrix of frequencies, one per each "cell"

ii = ceil((1:Size)/SizeCell);
FreqMat=Freq(ii,ii); %matrix of frequencies, all those within the same box are equal
%}

%% Generate the video

tic;
t = (1:N_frames)/FrameRate;
frame_stack = arrayfun(@(t) uint8(CIMat + AmplMat.* sin( 2*pi*FreqMat*t + PhaseMat)), t, 'UniformOutput',false );
frame_stack = cat(3,frame_stack{:});


%% Analysis

if matlabpool('size') == 0
    matlabpool open
end

for sbc = length(FOVsizevect):-1:1
    
    FOVsize = FOVsizevect(sbc);
    N_FOVs = (Size/FOVsize)^2;
    max_q = FOVsize/2;
    qvec = 2*pi/FOVsize * (1:max_q);
    
    
    disp(['DDM with SizeBox = ',num2str(FOVsize)]);
    tic
    chopped_frame_stack = mat2cell(frame_stack, FOVsize*ones(Size/FOVsize,1), FOVsize*ones(Size/FOVsize,1), N_frames); %Divide the video into SizeBox x SizeBox x N_frames videos
    
    Frequency = zeros(N_FOVs, max_q);
    Amplitude = zeros(N_FOVs, max_q);
    Damping   = zeros(N_FOVs, max_q);
    
    if sbc~=length(FOVsizevect)
        parfor i=1:N_FOVs
            [Iqtau{i},err_Iqtau{i}] = DDM_core(chopped_frame_stack{i});
            [Frequency(i,:), Damping(i,:), Amplitude(i,:)]=fit_Iqtau_FluctMat_new(Iqtau{i},err_Iqtau{i});
        end
    else
        for i=1:N_FOVs
            [Iqtau{i},err_Iqtau{i}] = DDM_core(chopped_frame_stack{i});
            [Frequency(i,:), Damping(i,:), Amplitude(i,:)]=fit_Iqtau_FluctMat_new(Iqtau{i},err_Iqtau{i});
        end
    end
    
    results(sbc).Frequency = Frequency;
    results(sbc).Damping = Damping;
    results(sbc).Amplitude = Amplitude;
    results(sbc).FOVsize = FOVsize;
    results(sbc).Iqtau = Iqtau;
    results(sbc).err_Iqtau = err_Iqtau;
    results(sbc).qvec = qvec;
    
    clear Iqtau err_Iqtau
    
    toc
    
end

%% Plots

ColorSet = varycolor(numel(FOVsizevect));

hf = figure;

for i=1:numel(FOVsizevect)
    
    if ~isvector(results(i).Frequency)
        Frequency4plot = mean(results(i).Frequency);
        Amplitude4plot = mean(results(i).Amplitude);
        Damping4plot   = mean(results(i).Damping);
    else
        Frequency4plot = results(i).Frequency;
        Amplitude4plot = results(i).Amplitude;
        Damping4plot   = results(i).Damping;
    end
    
    if isfield(results,'qvec')
        qvec = results(i).qvec;
    else
        qvec = 2*pi/results(i).FOVsize * (1:results(i).FOVsize/2);
    end
    
    hsp1 = subplot(221); hold on
    plot(qvec,Frequency4plot,'Color',ColorSet(i,:));
    
    hsp2 = subplot(222); hold on
    plot(qvec,Amplitude4plot,'Color',ColorSet(i,:));
    
    hsp3 = subplot(223); hold on
    plot(qvec,Damping4plot,'Color',ColorSet(i,:));
%     errorbar(qvec,Damping4plot,std(results(i).Damping)./sqrt(size(results(i).Damping,1)),'Color',ColorSet(i,:));
    
    leg{i} = ['FOV Size = ',num2str(results(i).FOVsize)];
    
    D(i) = mean( Damping4plot  );
    err_D(i) = std(Damping4plot);
    
end

hsp4 = subplot(224); hold on
errorbar(FOVsizevect,D,err_D,'.r');


set(hsp1, 'XLim', [0, pi]);
xlabel(hsp1, 'q [px^{-1}]', 'FontSize', 14, 'FontWeight', 'b');
ylabel(hsp1, 'Frequency [Hz]', 'FontSize', 14, 'FontWeight', 'b');


set(hsp2, 'XLim', [0, pi], 'YScale', 'log');
xlabel(hsp2, 'q [px^{-1}]', 'FontSize', 14, 'FontWeight', 'b');
ylabel(hsp2, 'Amplitude [a.u.]', 'FontSize', 14, 'FontWeight', 'b');
l = legend(hsp2, leg);
set(l, 'Location', 'NorthEast');

set(hsp3, 'XLim', [0, pi]);
xlabel(hsp3, 'q [px^{-1}]', 'FontSize', 14, 'FontWeight', 'b');
ylabel(hsp3, '1/\tau_{\sigma} [s^{-1}]', 'FontSize', 14, 'FontWeight', 'b');

set(hsp4, 'XScale','log','YScale','log');
xlabel(hsp4, 'Size of the DDM window', 'FontSize', 14, 'FontWeight', 'b');
ylabel(hsp4, '1/\tau_{\sigma} @q~1.2px^{-1} [s^{-1}]', 'FontSize', 14, 'FontWeight', 'b');
%}
%%
% save('hexagonal_patches_smaller.mat')
% save('square_1pxpatches_smaller.mat')
save('new_hexagonal_patches_2.mat')
% save('new_square_patches.mat')

