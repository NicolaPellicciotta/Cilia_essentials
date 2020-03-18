%%

% fit_Iqtau_colloids.m fits exponentials to find the characteristic times,
% tau_q, and then 'calls' refit_Iqtau_colloids.m to produce the choice
% scroll, choice, and end plots

%%
%clear all; close all;

%load results.mat;
% figure(101);
% imagesc(Iqtau(10:end,:));
% xlabel('lag time','FontSize',14);
% ylabel('q','FontSize',14);
% colorbar
% shg
% pause


%% px2mum: multiply measure in px by px2mum to have the measure in mum

px2mum = 0.5263; %Ximeas 10x
% px2mum = 0.2644; %Ximeas 20x
% px2mum = 0.1319; %Ximeas 40x
% px2mum = 0.5816; %Grasshopper 10x

%%

%max_q = min(128,size(Iqtau,1));
% max_tau = min(256,size(Iqtau,2));
% max_tau = 512;

frame_size = 2*size(Iqtau,1);
max_q = size(Iqtau,1);      %trying to extract data from the whole Iqtau matrix. Tighten these 2 parameters if the fit gives errors.
max_tau = size(Iqtau,2);


if ~exist('frames2s','var')
    framerate = input('Enter the framerate of the video: ');
    frames2s = 1/framerate;
    px2mum = input('How many microns per pixel? ');
end
tau = frames2s*[1:max_tau]'; % lag times in seconds

% initialisation of variables for the fit
A = zeros(max_q,1);
tau_q = zeros(max_q,1);
B = zeros(max_q,1);
err_a = A;
err_tau_q = tau_q;
err_B = B;

% chi2_colloids = zeros(max_q,1);

% function to fit
ft = fittype('a*(1-exp(-tau/b))+c','independent','tau');

% definition of the array of q
q = 2*pi/frame_size*[1:max_q]'; %in px^-1

% initialisation structure array for fit results
% fits_array = [];
Iqtau = 1e4 .* Iqtau; %or the fit won't work properly
% figure(102);
for q_c=1:max_q     %q_c == counter on q
    
    y =  Iqtau(q_c,1:max_tau)'; %just extracting a row (function of tau) of the Iqtau matrix
%     y = 1e4 * Iqtau(q_c,1:max_tau)'; %just extracting a row (function of tau) of the Iqtau matrix
%     err_y = 1e4 * err_Iqtau(q_c,1:max_tau)';
    fo = fitoptions('Method','NonLinearLeastSquares','Lower',[0, 0, 0],'Upper',[Inf, Inf, y(1)],'StartPoint',[max(Iqtau(q_c,:)'), 1, 0]);%,'Weights',1./err_y.^2);   % options for the fit
    fit_out = fit(tau,y,ft,fo);
    A(q_c) = fit_out.a;             % saturation value
    tau_q(q_c) = fit_out.b;      % time constant (this is what this fit is done for)
    B(q_c) = fit_out.c;             % offset
    temp = 0.5*diff(confint(fit_out,0.682));    % this and next three lines: errors on A, tau_q and B
    err_A(q_c) = temp(1);
    err_tau_q(q_c) = temp(2);
    err_B(q_c) = temp(3);
    
%     chi2_colloids(q_c) = sum( (y - A(q_c)*(1-exp(-tau/tau_q(q_c)))+B(q_c)).^2 ./ err_y.^2 ) / (max_tau-4);
    
%     plot stuff
%     cla;
%     plot(tau,y,'.');
%     hold on;
%     plot(tau, A(q_c)*(1-exp(-tau/tau_q(q_c)))+B(q_c), 'r');
%     plot(tau, fit_out(tau), 'g');
%     title(['mode = ',num2str(q_c),',  \lambda = ',num2str(frame_size/q_c),'px'],'FontSize',14);
%     xlabel('Lag time \tau [s]', 'FontSize',14);
%     ylabel('Intensity [a.u.]', 'FontSize',14);
%     legend('experimental points','fit')
%     shg
%     pause(0.01);
%     

    fits_array{q_c} = fit_out;
    
end

Iqtau = 1e-4 .* Iqtau; % prevent Iqtau tending to infinity


%%

refit_Iqtau_colloids

%%
% scroll_through_Iqtau_fits(Iqtau,fits_array,tau);
% 
% %% fit of the time constant tau_q vs q
% 
% % change the next two lines if you want to fit a different range (in q) of the datapoints
% min_q_fit = 1;
% max_q_fit = length(q);
% 
% % plot stuff
% clf;
% hold on;
% % errorbar(q(min_q_fit:max_q_fit),tau_q(min_q_fit:max_q_fit),err_tau_q(min_q_fit:max_q_fit),'g.');
% % dummy_q = [min_q_fit:max_q_fit]';                                                     %use this line instead of the following if you get an unexpected error
% dummy_q = [1:length(q)]';
% errorbar(dummy_q(min_q_fit:max_q_fit),tau_q(min_q_fit:max_q_fit),err_tau_q(min_q_fit:max_q_fit),'g.');
% set(gca,'XScale','log','YScale','log');
% % xlabel('q [px^{-1}]','FontSize',14);
% xlabel('mode','FontSize',14);
% ylabel('\tau_q [s]','FontSize',14);
% 
% % interactive stuff
% disp('Click on the leftmost and on the rightmost of the points you want to fit, then press enter');
% title('Click on the leftmost and on the rightmost of the points you want to fit, then press enter');
% [x,y] = getpts(gcf);
% if x(2)<x(1)
%     x = x(end:-1:1);
%     y = y(end:-1:1);
% end
% % [~,min_q_fit]=min((x(1)-q).^2+(y(1)-tau_q).^2);
% % [~,max_q_fit]=min((x(2)-q).^2+(y(2)-tau_q).^2);
% [~,min_q_fit]=min((x(1)-dummy_q).^2+(y(1)-tau_q).^2);
% [~,max_q_fit]=min((x(2)-dummy_q).^2+(y(2)-tau_q).^2);
% %%
% % plot stuff
% figure(103);
% clf;
% hold on;
% errorbar(q([1:min_q_fit,max_q_fit:end]),tau_q([1:min_q_fit,max_q_fit:end]),err_tau_q([1:min_q_fit,max_q_fit:end]),'b.');
% errorbar(q(min_q_fit:max_q_fit),tau_q(min_q_fit:max_q_fit),err_tau_q(min_q_fit:max_q_fit),'g.');
% set(gca,'XScale','log','YScale','log');
% xlabel('q [px^{-1}]','FontSize',14);
% ylabel('\tau_q [s]','FontSize',14);
% 
% 
% % actual fit
% foo = fit(q(min_q_fit:max_q_fit),tau_q(min_q_fit:max_q_fit),fittype('1/(Dm*x^2)','independent','x'),fitoptions('Method','NonLinearLeastSquares','Weights',1./err_tau_q(min_q_fit:max_q_fit).^2));
% foo2 = fit(q(min_q_fit:max_q_fit),tau_q(min_q_fit:max_q_fit),fittype('1/(Dm*x^alpha)','independent','x'),fitoptions('Method','NonLinearLeastSquares','Weights',1./err_tau_q(min_q_fit:max_q_fit).^2));
% 
% % plot stuff
% plot(q(min_q_fit:max_q_fit),1./(foo.Dm*q(min_q_fit:max_q_fit).^2),'r');
% plot(q(min_q_fit:max_q_fit),1./(foo2.Dm*q(min_q_fit:max_q_fit).^foo2.alpha),'k');
% set(gca,'XScale','log','YScale','log');
% xlabel('q [px^{-1}]','FontSize',14);
% ylabel('\tau_q [s]','FontSize',14);
% legend('non-fitted data points','fitted data points',['fit with free exponent \alpha = ',num2str(foo2.alpha)],'fit with exponent fixed at 2','Location','SouthWest')
% % conversion to mum2/s
% Dm = foo.Dm * px2mum^2;       % now in mum2/s (foo.Dm is in px2/s)
% DmFo = foo2.Dm * px2mum^2;
% disp(['Diffusion Coefficient (Fitted) = ', num2str(Dm),' um2/s']);
% disp(['Diffusion Coefficient (Forced) = ', num2str(DmFo),' um2/s']);



% Iqtau = 1e-4 .* Iqtau; % prevent Iqtau tending to infinity