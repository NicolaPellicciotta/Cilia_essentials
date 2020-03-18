warning('This will take a very long time and overwrite saved data')
decision = input('Are you sure you want to run the whole script? ');

if ~decision,
    return
end

disp('blabla')

clear all
close all

%% structure initialisation

CommonData.Size = 480; %in pixel, size of the video
CommonData.SizeBox = 1; %size of the box of pixels with the same frequency. Now 1 = 1px
CommonData.N_box = CommonData.Size/CommonData.SizeBox; % Number of boxes of pixels with the same frequency
CommonData.TotalTime = 15; % in seconds
CommonData.FrameRate = 300; %fps
CommonData.N_frames =  CommonData.TotalTime * CommonData.FrameRate;
CommonData.Time = linspace(1/CommonData.FrameRate, CommonData.TotalTime, CommonData.N_frames);
CommonData.max_n = CommonData.Size / 2;
CommonData.max_tau = CommonData.N_frames / 2;

%% data structure initialisation

sigvec = linspace(0,0.3,31); %values of sigma that will be investigated
pvec = 0.0:0.05:0.5; %values of probabilities of phase-slip that will be investigated

ii = ceil( (1 : CommonData.Size) / CommonData.SizeBox); %index to fill the FreqMat giving to each pixel the frequency of the right box. only useful when SizeBox ~= 1;

for p = numel(pvec):-1:1;
    for s = numel(sigvec):-1:1;
        SimData(s,p).f = 15; %Hertz
        SimData(s,p).sigma = sigvec(s);
        SimData(s,p).pPS = pvec(p);
        SimData(s,p).Freq = normrnd(SimData(s,p).f, SimData(s,p).sigma, CommonData.N_box .* [1 1]); % matrix, one frequency per each box
        SimData(s,p).FreqMat = SimData(s,p).Freq(ii,ii); %matrix of frequencies, all those within the same box are equal, but is Size x Size
        SimData(s,p).PhaseMat = 2 * pi * rand(CommonData.N_box);
        SimData(s,p).PhaseSlipAmp = pi/32;
        SimData(s,p).Iqtau = zeros(CommonData.max_n, CommonData.max_tau);
        SimData(s,p).err_Iqtau = zeros(CommonData.max_n, CommonData.max_tau);
        SimData(s,p).Envelope = [];
        SimData(s,p).TimesEnvelope = [];
        SimData(s,p).StdEnvelope = [];
        SimData(s,p).NormEnvelope = [];
        SimData(s,p).NormStdEnvelope = [];
        SimData(s,p).FitNormEnvelope = [];
        SimData(s,p).AmS = NaN;
        SimData(s,p).DTS = NaN;
        SimData(s,p).AmP = NaN;
        SimData(s,p).DTP = NaN;
        SimData(s,p).Off = NaN;
    end
end

%% "video" creation

AmplMat = 127*ones(CommonData.Size); %for uint8 images
CIMat = 128*ones(CommonData.Size);

pool = gcp; %open pool of workers

control = false(numel(SimData),1);

for ccc = 1:6
    
    parfor cc = (ccc-1) * 60 + 1 : min( (ccc) * 60 + 1, numel(SimData)) % Configuration Counter
        
        control(cc) = true; %control to have swept all the parameters
        
        % initialisation of variables for video creation
        PhaseMat = SimData(cc).PhaseMat;
        frame_stack = zeros(CommonData.Size,CommonData.Size,CommonData.N_frames,'uint8');
        
        % here starts the actual video creation
        for tt = 1 : CommonData.N_frames %frame counter
            t = tt/CommonData.FrameRate;   %real time, in seconds
            p = rand(CommonData.Size) < SimData(cc).pPS; %pixels that undergo a phase-slip
            PhaseMat(p) = PhaseMat(p) + (1-2*round(rand(sum(p(:)),1))) * SimData(cc).PhaseSlipAmp; % increase/decrease of the phase
            frame_stack(:,:,tt) = uint8(CIMat + AmplMat .* sin( 2 * pi * SimData(cc).FreqMat * t + PhaseMat  ));
        end
        
        % DDM analysis
        [SimData(cc).Iqtau, SimData(cc).err_Iqtau] = DDM_core(frame_stack);
        [SimData(cc).Envelope, locs] = findpeaks( mean(SimData(cc).Iqtau), 'minpeakdistance', round( CommonData.FrameRate / (0.7 * SimData(cc).f) ));
        SimData(cc).TimesEnvelope = CommonData.Time(locs);
        dummy = std(SimData(cc).Iqtau);
        SimData(cc).StdEnvelope = dummy(locs)./sqrt(CommonData.max_n);
    end
    
    
end





%% Normalisation

minvalue = min(horzcat(SimData(:).Envelope));
for cc = 1:numel(SimData)
    SimData(cc).NormEnvelope = SimData(cc).Envelope - minvalue;
end
for cc = 1:numel(SimData)
    maxvalue = max(SimData(cc).NormEnvelope);
    SimData(cc).NormEnvelope = SimData(cc).NormEnvelope ./ maxvalue;
    SimData(cc).NormStdEnvelope = SimData(cc).StdEnvelope ./ maxvalue;
end

figure(202);
hold on;
xlim([-1 8]);
ylim([-eps 1+eps]);

for s = 1:numel(sigvec)
    for p = 1:numel(pvec)
        if SimData(s,p).sigma ~= 0 || SimData(s,p).pPS ~= 0
            
            tt = SimData(s,p).TimesEnvelope;
            yy = SimData(s,p).NormEnvelope;
            eyy = SimData(s,p).NormStdEnvelope;
            
            
            ft = fittype('a * exp( -b * x - c^2 * x^2 ) + off');
            
            fo = fitoptions('Method','NonLinearleastSquares',...
                'Lower', [0 0 0 0], ...
                'Upper', [1 Inf Inf 0.01],...
                'StartPoint',[0.9  0.8 0.8 0],...
                'MaxFunEvals', 1e7, ...
                'MaxIter', 1e7, ...
                'TolFun', 1e-14,...
                'Weights',1./eyy.^2);
            
            fit_out = fit(tt',yy',ft,fo);
            SimData(s,p).AmS = fit_out.a;
            SimData(s,p).DTS = 1./fit_out.c;
            SimData(s,p).DTP = 1./fit_out.b;
            SimData(s,p).Off = fit_out.off;
            
        end
        
        cla
        errorbar(SimData(s,p).TimesEnvelope, SimData(s,p).NormEnvelope, SimData(s,p).NormStdEnvelope ,'.');
        if SimData(s,p).sigma ~= 0 || SimData(s,p).pPS ~= 0
            SimData(s,p).FitNormEnvelope = fit_out(SimData(s,p).TimesEnvelope)';
            plot(SimData(s,p).TimesEnvelope, SimData(s,p).FitNormEnvelope,'r');
        end
        title(['\sigma = ',num2str(SimData(s,p).sigma),'    pPS = ',num2str(SimData(s,p).pPS)]);
        pause(0.01);
        
        
    end
end



%% uberplot log envelope vs tau2
figure(203);

set(gcf,'Position',get(0,'ScreenSize')*8/10,'Units','Normalized');
margin = 0.07;

S = 3;
P = 3;

hs = zeros(P,S);
hp = zeros(P,S);
he = zeros(P,S);

for s = S:-1:1
    for p = P:-1:1
        
        hs(p,s) = subplot('Position',[ margin + (1-2*margin)/S * (s-1) , margin + (1-2*margin)/P * (p-1), (1-2*margin)/S, (1-2*margin)/P],...
            'Box','on',...
            'XLim',[-0.5 7.5].^2,...
            'YLim',[-4.5 0.8],...
            'XTickLabel',[],...
            'YTickLabel',[]);
        
        if p == P
            set(gca,'XAxisLocation','top');
            xlabel(num2str(SimData(2*s-1,2*p-1).sigma),'FontSize',18);
        end
        
        if s == 1
            set(gca,'YAxisLocation','left');
            ylabel(num2str(SimData(2*s-1,2*p-1).pPS),'FontSize',18);
        end
        
        if p == (P+1)/2 && s ==S
            set(gca,'YAxisLocation','right');
            ylabel('Phase-slip probability','FontSize',18,'FontWeight','b');
        end
        
        if p == 1
            set(gca,'XAxisLocation','bottom');
            if s == (S+1)/2
                xlabel('\sigma of frequency distribution, [Hz]','FontSize',18,'FontWeight','b');
            end
            if s == S
                xlabel('\tau^2, [s^2]','FontSize',16);
                ylabel('Amplitude, [a.u.]','FontSize',16);
                set(gca,'XTickLabel',get(gca,'XTick'));
                set(gca,'YAxisLocation','right');
                set(gca,'YTickLabel',sprintf('%.2f|',exp(get(gca,'YTick'))));               
            end
        end
                
%         her(p,s) = errorbar(SimData(2*s-1,2*p-1).TimesEnvelope.^2, log(SimData(2*s-1,2*p-1).NormEnvelope), SimData(2*s-1,2*p-1).NormStdEnvelope ./ SimData(2*s-1,2*p-1).NormEnvelope ,'.','MarkerSize',10);
        hold on;
        he(p,s) = plot(SimData(2*s-1,2*p-1).TimesEnvelope.^2, log(SimData(2*s-1,2*p-1).NormEnvelope),'.','MarkerSize',10);
        if SimData(2*s-1,2*p-1).sigma ~= 0 || SimData(2*s-1,2*p-1).pPS ~= 0
            hp(p,s) = plot(SimData(2*s-1,2*p-1).TimesEnvelope.^2, log(SimData(2*s-1,2*p-1).FitNormEnvelope),'r','LineWidth',1.2);
        end
        
        if p==P && s ==S
            leg = legend('envelope of <I(q,\tau)>_q','fit');
            set(leg,'Location','NorthEast','box','off','FontSize',14);
        end
        
        set(hs(p,s),'XLim',[-0.5 7.5].^2,...
            'YLim',[-4.5 0.8]);
        
    end
end

set(gcf,'Units','Normalized');
set(gcf,'PaperPositionMode','auto');
set(gcf,'InvertHardcopy','off','Color','w');
saveas(gcf,'logEnvelope_vs_tau2.eps','epsc2')
export_fig 'logEnvelope_vs_tau2' -png -a1 -native
export_fig 'logEnvelope_vs_tau2' -pdf -q101
plot2svg('logEnvelope_vs_tau2.svg',gcf)



%% uberplot log envelope vs tau

figure(204); % log envelope vs tau

set(gcf,'Position',get(0,'ScreenSize')*8/10,'Units','Normalized');
margin = 0.07;

S = 3;
P = 3;

hs = zeros(P,S);
hp = zeros(P,S);
he = zeros(P,S);

for s = S:-1:1
    for p = P:-1:1
        
        hs(p,s) = subplot('Position',[ margin + (1-2*margin)/S * (s-1) , margin + (1-2*margin)/P * (p-1), (1-2*margin)/S, (1-2*margin)/P],...
            'Box','on',...
            'XLim',[-0.5 7.5],...
            'YLim',[-4.5 0.8],...
            'XTickLabel',[],...
            'YTickLabel',[]);
        
        if p == P
            set(gca,'XAxisLocation','top');
            xlabel(num2str(SimData(2*s-1,2*p-1).sigma),'FontSize',18);
        end
        
        if s == 1
            set(gca,'YAxisLocation','left');
            ylabel(num2str(SimData(2*s-1,2*p-1).pPS),'FontSize',18);
        end
        
        if p == (P+1)/2 && s ==S
            set(gca,'YAxisLocation','right');
            ylabel('Phase-slip probability','FontSize',18,'FontWeight','b');
        end
        
        if p == 1
            set(gca,'XAxisLocation','bottom');
            if s == (S+1)/2
                xlabel('\sigma of frequency distribution, [Hz]','FontSize',18,'FontWeight','b');
            end
            if s == S
                xlabel('\tau, [s]','FontSize',16);
                ylabel('Amplitude, [a.u.]','FontSize',16);
                set(gca,'XTickLabel',get(gca,'XTick'));
                set(gca,'YAxisLocation','right');
                set(gca,'YTickLabel',sprintf('%.2f|',exp(get(gca,'YTick'))));               
            end
        end
                
        %         he(p,s) = errorbar(SimData(s,p).TimesEnvelope, SimData(s,p).NormEnvelope, SimData(s,p).NormStdEnvelope ,'.');
        hold on;
        he(p,s) = plot(SimData(2*s-1,2*p-1).TimesEnvelope, log(SimData(2*s-1,2*p-1).NormEnvelope),'.','MarkerSize',10);
        if SimData(2*s-1,2*p-1).sigma ~= 0 || SimData(2*s-1,2*p-1).pPS ~= 0
            hp(p,s) = plot(SimData(2*s-1,2*p-1).TimesEnvelope, log(SimData(2*s-1,2*p-1).FitNormEnvelope),'r','LineWidth',1.2);
        end
        
        if p==P && s ==S
            leg = legend('envelope of <I(q,\tau)>_q','fit');
            set(leg,'Location','NorthEast','box','off','FontSize',14);
        end
    end
end

set(gcf,'Units','Normalized');
set(gcf,'PaperPositionMode','auto'); %qui penso ci sia da dargli la misura in cm
set(gcf,'InvertHardcopy','off','Color','w');
saveas(gcf,'logEnvelope_vs_tau.eps','epsc2')
export_fig 'logEnvelope_vs_tau' -png -a1 -native
export_fig 'logEnvelope_vs_tau' -pdf -q101
plot2svg('logEnvelope_vs_tau.svg',gcf)



%% fit results map => plot map_decay_times

AmS = reshape([SimData(:).AmS],size(SimData));
DTS = reshape([SimData(:).DTS],size(SimData));
DTP = reshape([SimData(:).DTP],size(SimData));
Off = reshape([SimData(:).Off],size(SimData));
ss = reshape([SimData(:).sigma],size(SimData));
pp = reshape([SimData(:).pPS],size(SimData));

figure(205);

% subplot(2,1,1)
% 
% hax = gca;
% set(hax,'YDir','reverse');
% hams = imagesc(rot90(AmS));
% set(hax,'XTick',1:5:numel(sigvec));
% set(hax,'XTickLabel',sigvec(get(hax,'XTick')));
% set(hax,'YTick',1:2:numel(pvec));
% set(hax,'YTickLabel',pvec(fliplr(get(hax,'YTick'))));
% axis image;

subplot(2,1,1)

hax = gca;
set(hax,'YDir','reverse');
hdts = imagesc(1./rot90(DTS));
set(hax,'XTick',1:5:numel(sigvec));
set(hax,'XTickLabel',sigvec(get(hax,'XTick')));
set(hax,'YTick',1:2:numel(pvec));
set(hax,'YTickLabel',pvec(fliplr(get(hax,'YTick'))));
axis image;
xlabel('\sigma of frequency distribution, [Hz]','FontSize',14,'FontWeight','b');
ylabel({'Phase-slip', 'probability'},'FontSize',14,'FontWeight','b');
hc = colorbar;
ylabel(hc,'1/\tau_{\sigma}, [s^{-1}]','FontSize',14,'FontWeight','b');


subplot(2,1,2)

hax = gca;
set(hax,'YDir','reverse');
hdtp = imagesc(1./(rot90(DTP)));
set(hax,'XTick',1:5:numel(sigvec));
set(hax,'XTickLabel',sigvec(get(hax,'XTick')));
set(hax,'YTick',1:2:numel(pvec));
set(hax,'YTickLabel',pvec(fliplr(get(hax,'YTick'))));
axis image;
xlabel('\sigma of frequency distribution, [Hz]','FontSize',14,'FontWeight','b');
ylabel({'Phase-slip', 'probability'},'FontSize',14,'FontWeight','b');
hc = colorbar;
ylabel(hc,'1/\tau_p, [s^{-1}]','FontSize',14,'FontWeight','b');

shg


set(gcf,'Units','Normalized');
set(gcf,'PaperPositionMode','auto');
set(gcf,'InvertHardcopy','off','Color','w');
saveas(gcf,'map_decay_times.eps','epsc2')
export_fig 'map_decay_times' -png -a1 -native
export_fig 'map_decay_times' -pdf -q101

%% plot fit decay_time vs pPS and sigma

figure(206), clf
centred_resized_figure(gcf, [1 2]);

subplot(121)
hold on
% hpp = errorbar(pvec, nanmean(1./DTP), nanstd(1./DTP),'o','LineWidth',1.5);
hpp = plot(pvec, nanmean(1./DTP),'o','MarkerSize',6,'LineWidth',1.2);
polyfit(pvec, nanmean(1./DTP),1)
hfp = plot(pvec,polyval(polyfit(pvec, nanmean(1./DTP), 1),pvec),'r','LineWidth',1.2);
xlabel('Phase-slip probability','FontSize',14,'FontWeight','b');
ylabel('1/\tau_p, [s^{-1}]','FontSize',14,'FontWeight','b');
xlim([pvec(1)-0.01, pvec(end)+0.01])
ylim([- 0.01, 0.701])
box on
hlp = legend('data','linear fit','Location','NW');
set(hlp,'Box','off');

subplot(122)
hold on
% hps = errorbar(sigvec,nanmean(1./DTS,2), nanstd(1./DTS,[],2),'o','LineWidth',1.5);
hps = plot(sigvec,nanmean(1./DTS,2),'o','MarkerSize',6,'LineWidth',1.2);
polyfit(sigvec', nanmean(1./DTS,2),1)
hfs = plot(sigvec,polyval(polyfit(sigvec', nanmean(1./DTS,2),1),sigvec),'r','LineWidth',1.2);
xlabel('\sigma of frequency distribution, [Hz]','FontSize',14,'FontWeight','b');
ylabel('1/\tau_{\sigma}, [s^{-1}]','FontSize',14,'FontWeight','b');
xlim([sigvec(1)-0.007, sigvec(end)+0.01])
ylim([- 0.01, 1.401])
box on
hls = legend('data','linear fit','Location','NW');
set(hls,'Box','off');

set(gcf,'Units','Normalized');
set(gcf,'PaperPositionMode','auto');
set(gcf,'InvertHardcopy','off','Color','w');
% saveas(gcf,'av_decay_times_vs_sigma_pPS.eps','epsc2')
% export_fig 'av_decay_times_vs_sigma_pPS' -png -a1 -native
% export_fig 'av_decay_times_vs_sigma_pPS' -pdf -q101
% plot2svg('av_decay_times_vs_sigma_pPS.svg',gcf)


 


%% save

save parsweep_sigma_pPS_long_try -v7.3;


















