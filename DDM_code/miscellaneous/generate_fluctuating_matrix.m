clear all
close all

% sigvec = logspace(log10(0.1),log10(0.2),10);
pvec = 1-logspace(log10(0.01),log10(0.10),30);
pvec = 1-0.4;
for pc = 1:length(pvec)
    
    %% Constants
    Size = 480;
    SizeBox = 1;
    TotalTime = 10; % in seconds
    FrameRate = 300; %fps
    
    f = 15;  %mean frequency
    sigma = 0; %sigma
    
    pPS = 1-pvec(pc)
    
    %% Parameters initialisation
    AmplMat = 128*ones(Size); %for uint8 images
    CIMat = 128*ones(Size);
    % PhaseMat = repmat(linspace(0,2*pi,Size),[Size,1]);
    PhaseMat = zeros(Size);
    PhaseNoiseAmp = 0;
    
    PhaseSlipAmp = pi/32;
    
    Freq = f+sigma*randn(Size/SizeBox); %matrix of frequencies, one per each box
    
    ii = ceil((1:Size)/SizeBox);
    FreqMat=Freq(ii,ii); %matrix of frequencies, all those within the same box are equal
    
    %% Time evolution
    
    
    frame_stack = zeros(Size,Size,TotalTime,'uint8');
    
    for tt=1:TotalTime*FrameRate %time counter, in frames
        
        t = tt/FrameRate;   %real time, in seconds
        p = rand(Size)>pvec(pc); %probability of phase slip
        PhaseMat(p) = PhaseMat(p) + (1-2*round(rand(sum(p(:)),1))) * PhaseSlipAmp;
        frame_stack(:,:,tt) = uint8(CIMat + AmplMat.* sin( 2*pi*FreqMat*t + PhaseMat + PhaseNoiseAmp*(0.5-rand) ));
        
    end
    
    [Iqtau, err_Iqtau] = DDM_core(frame_stack);
    
    save(['FluctMatPhaseSlipPix32e-1_Size=',num2str(Size),'_Box=',num2str(SizeBox),...
        '_freq=',num2str(f),'_sigma=',num2str(sigma),'_RunTime=',...
        num2str(TotalTime),'_fps=',num2str(FrameRate),'_pPS=',num2str(pPS),'.mat']);
    
    clearvars -except pvec pc
    
    
    
end

beep