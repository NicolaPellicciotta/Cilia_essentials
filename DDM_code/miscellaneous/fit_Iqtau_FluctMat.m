function [ Frequency, Damping, Amplitude, gof] = fit_Iqtau_FluctMat( Iqtau, err_Iqtau, max_q, max_tau )
%fit_Iqtau_FluctMat Fits every row of the Iqtau matrix with a dampened
%oscillatory function. works on synthetic data
%   Detailed explanation goes here

if nargin < 4 || isempty(max_tau)
    max_tau = size(Iqtau,2);
end
if nargin < 3 || isempty(max_q)
    max_q = size(Iqtau,1);
end


Amplitude = zeros(max_q, 1, 'double');
Frequency = zeros(max_q, 1, 'double');
Damping   = zeros(max_q, 1, 'double');
Offset    = zeros(max_q, 1, 'double');

err_Amplitude = zeros(max_q, 1, 'double');
err_Frequency = zeros(max_q, 1, 'double');
err_Damping   = zeros(max_q, 1, 'double');
err_Offset    = zeros(max_q, 1, 'double');

FrameRate = 300; %fps

for qq = 1:max_q
    
    xx = (1:max_tau)'/FrameRate;
    yy = Iqtau(qq,1:max_tau)';
%     ww = 1./err_Iqtau(qq,1:max_tau)'.^2;
    ww = 1./xx;
    
    %% initialisation of the fit parameters
    
    init_Offset = mean(yy);
    init_Amplitude = max(yy)+init_Offset;
    init_Damping = 1e5;
    init_Frequency = 15; %Hertz, known
    
    fit_options = fitoptions('Method','NonLinearleastSquares',...
        'Startpoint', [init_Amplitude  init_Damping init_Frequency init_Offset], ...
        'Lower', [0 0 0 0], ...
        'Upper', [(max(yy)-min(yy)) Inf Inf max(yy) ], ...
        'MaxFunEvals', 1e4, ...
        'MaxIter', 1e4, ...
        'TolFun', 1e-5*max(yy), ...
        'Weights',ww);

    fit_type = fittype('Ampl * cos(2*pi*Freq*x + pi) * exp( - ( 2*pi*Freq*x )^2 / Damp)+Offset');
    
    [fit_output, gof, fit_info] = fit(xx, yy, fit_type, fit_options);
    
       
    
    Amplitude(qq) = fit_output.Ampl;
    Frequency(qq) = fit_output.Freq;
    Damping(qq)   = fit_output.Damp;
    Offset(qq)    = fit_output.Offset;
    
    %         errors = 0.5*diff(confint(fit_output));
    %
    %         err_Amplitude(qq) = errors(1);
    %         err_Frequency(qq) = errors(2);
    %         err_Damping(qq)   = errors(3);
    %         err_Offset(qq)    = errors(4);
    
    
    
    clf;
    plot(xx,yy,'.');
    hold on
    plot(xx, fit_output(xx),'r');
    title(['q = ',num2str(qq)])
    pause(0.01)
    
    
end

end

