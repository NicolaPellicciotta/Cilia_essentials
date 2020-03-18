function [cilia_fits,R] = fit_Iqtau_cilia_nicola( Iqtau, max_q, max_tau )
%fit_Iqtau_cilia Fits every row of the Iqtau matrix with a dampened
%oscillatory function
%   Detailed explanation goes here

if nargin < 3 || isempty(max_tau)
    max_tau = size(Iqtau,2);
end
if nargin < 2 || isempty(max_q)
    max_q = size(Iqtau,1);
end


Amplitude = nan(max_q, 1, 'double');
Frequency = nan(max_q, 1, 'double');
Damping   = nan(max_q, 1, 'double');
Offset    = nan(max_q, 1, 'double');
GOF       = nan(max_q, 1, 'double');

%% preliminary search for qstart
%|    _
% \  / \
%  \/   '.__
% Since the shape of the <Iqtau>tau is as above, looking for the negative
% peak and then mediating over q only from then on
% 
% [~,qstart]=findpeaks(-1*mean(Iqtau(1:max_q,1:max_tau),2),'NPeaks',1);
% if isempty(qstart) || qstart > max_q/2
%     qstart = 1;
% end

qstart=1;
%% preliminary fit of the <Iqtau>q to get good parameters for the fit of each I(q',tau)
xx = (1:max_tau)';
yy = mean(Iqtau(qstart:max_q,1:max_tau))';

init_Offset = mean(yy(round(end/2):end));
init_Amplitude = abs(-yy(1)+init_Offset);
init_DecayTime = length(yy)/8;

autocorr = xcorr(yy-init_Offset);
dummy_freq = fftshift(fft( autocorr(length(yy):end) ));
dummy_freq = dummy_freq(ceil(length(dummy_freq)/2)+1:end);
[pk, init_f]=findpeaks(abs(dummy_freq));
init_f = init_f(pk == max(pk));
if isempty(init_f);init_f=1.1;end
    
init_Freq = 2*pi*init_f/length(autocorr(length(yy):end));





ft = fittype('Ampl*(1-exp(cos(Freq*xx)))*exp(-Damp*xx)+Offset','Independent','xx');
fo = fitoptions('Method','NonLinearLeastSquare',...
    'StartPoint',[init_Amplitude, 1/init_DecayTime, init_Freq, init_Offset],...
    'Lower',[0 0 abs(init_Freq)/2 0],...
    'Upper',[max(yy)+init_Offset, Inf, abs(init_Freq)*2, max(yy)],...
    'MaxFunEvals',1e6,'MaxIter',1e6,'TolFun',1e-9*init_Amplitude,'Weights',1./sqrt(xx));

alt_fo = fitoptions('Method','NonLinearLeastSquare',... %alternative fitoptions in case the first ones yield error
    'StartPoint',[init_Amplitude, 1/init_DecayTime, init_Freq, init_Offset],...
    'MaxFunEvals',1e6,'MaxIter',1e6,'TolFun',1e-9*init_Amplitude,'Weights',1./sqrt(xx));

try
    [fit_out] = fit(xx,yy,ft,fo);
catch
    findpeaks(-1*mean(Iqtau(1:max_q,1:max_tau),2),'NPeaks',1)
    [fit_out] = fit(xx,yy,ft,alt_fo);
end

init_Frequency = fit_out.Freq;
if init_Frequency <= 0, init_Frequency = init_Freq; end
init_Damping = fit_out.Damp;

%     figure
%     plot(xx,yy,'.');hold on
%     plot(xx,cilia_fits(qq).fit_out(xx),'r');
%     pause

%% fit of the I(q',tau) curves
%hpool = gcp; %if no parpool opened, opens a new one
% parfor qq = 1:max_q
for qq = 1:max_q
    
         
    % normalising y to make fit work better
    xx = (1:max_tau)';
    yy = Iqtau(qq,1:max_tau)';
    max_yy = max(yy);
    yy = yy./ max_yy;

    
    %% initialisation of the fit parameters
    
    init_Offset = mean(yy(round(end/2):end));
    init_Amplitude = abs(-yy(1)+init_Offset);
    %     init_Damping = length(yy)/4;
    
    
    %%% not enough signal if
    if init_Amplitude < 0.03*init_Offset, cilia_fits(qq).fit_out=nan;cilia_fits(qq).norm_fact=nan;end 
    
    if init_Amplitude < 0.03*init_Offset, continue; end;
%    disp(qq);
    fo = fitoptions('Method','NonLinearLeastSquare',...
        'StartPoint',[init_Amplitude, init_Damping, init_Frequency, init_Offset],...
        'Lower',[0 0 abs(init_Frequency)/2 0],...
        'Upper',[max(yy)+init_Offset, Inf, abs(init_Frequency)*2, max(yy)],...
        'MaxFunEvals',1e6,'MaxIter',1e6,'TolFun',1e-9*init_Amplitude,'Weights',1./sqrt(xx));
    
    alt_fo = fitoptions('Method','NonLinearLeastSquare',... %alternative fitoptions in case the first ones yield error
        'StartPoint',[init_Amplitude, init_Damping, init_Frequency, init_Offset],...
        'MaxFunEvals',1e6,'MaxIter',1e6,'TolFun',1e-9*init_Amplitude,'Weights',1./sqrt(xx));
    
    try
        [cilia_fits(qq).fit_out] = fit(xx,yy,ft,fo);
    catch
        [cilia_fits(qq).fit_out] = fit(xx,yy,ft,alt_fo);
    end
    
    
    cilia_fits(qq).norm_fact = max_yy;
    
    Amplitude(qq) = cilia_fits(qq).fit_out.Ampl;
    Frequency(qq) = abs(cilia_fits(qq).fit_out.Freq);
    Damping(qq)   = cilia_fits(qq).fit_out.Damp;
    Offset(qq)    = cilia_fits(qq).fit_out.Offset;
    GOF(qq) = sum(  abs(yy(1:floor(end/4))-cilia_fits(qq).fit_out(xx(1:floor(end/4)))) ) ;
    

    
    
end

Frequency = Frequency/(2*pi);
% 
R.Amplitude = Amplitude;
R.Frequency= Frequency;
R.Damping= Damping;
R.Offset = Offset;
R.GOF =GOF;
% 100 is the framerate used for Cedar videos


end

