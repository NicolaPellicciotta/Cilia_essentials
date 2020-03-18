% clear all; close all;

% filename = 'ballistic_colloids_v=50_s=0_TumbleProb=0.05.mat';

% load(filename);
% figure;
% imagesc(Iqtau(10:end,:));
% xlabel('lag time','FontSize',14);
% ylabel('q','FontSize',14);
% colorbar
% shg
% pause


%% px2mum: multiply measure in px by px2mum to have the measure in mum

% px2mum = 0.5263; %Ximeas 10x
% px2mum = 0.2644; %Ximeas 20x
% px2mum = 0.1319; %Ximeas 40x
% px2mum = 0.5816; %Grasshopper 10x

%%

%max_q = min(128,size(Iqtau,1));
% max_tau = min(256,size(Iqtau,2));
% max_tau = 512;

frame_size = 2*size(Iqtau,1);
max_q = size(Iqtau,1)/2;      %trying to extract data from the whole Iqtau matrix. Tighten these 2 parameters if the fit gives errors.
% max_q = 20;      %trying to extract data from the whole Iqtau matrix. Tighten these 2 parameters if the fit gives errors.
max_tau = size(Iqtau,2)/2;
% frames2s = 1/100;
% px2mum = 1;
if ~exist('frames2s','var')
    framerate = input('Enter the framerate of the video: ');
    frames2s = 1/framerate;
    px2mum = input('How many microns per pixel? ');
end
tau = frames2s*[1:max_tau]'; % lag times in seconds

% initialisation of variables for the fit
A = zeros(max_q,1);
qvel = zeros(max_q,1);
B = zeros(max_q,1);
err_a = A;
err_qvel = qvel;
err_B = B;

% chi2_ballistic = zeros(max_q,1);

% function to fit
ft = fittype('a*(1-sin(b*tau)/(b*tau))+c','independent','tau');

% definition of the array of q
q = 2*pi/frame_size*[1:max_q]'; %in px^-1

figure;
for q_c=1:max_q     %q_c == counter on q
    
    %     max_tau = size(Iqtau,2)-10*q_c;
    tau = frames2s*[1:max_tau]'; % lag times in seconds
    
    
    y = Iqtau(q_c,1:max_tau)'; %just extracting a row (function of tau) of the Iqtau matrix
    y=y*1e4;
%     err_y = 1e4* err_Iqtau(q_c,1:max_tau)';
    %     err_y = [0.02*max(y)*ones(1,round(max_tau/4)),0.2*max(y)*ones(1,max_tau-round(max_tau/4))]';
    
    dummy=abs(fft(y));
    dummy=dummy(2:round(max_tau/2));
    [~,wh]=max(dummy);
    start_b = 2*pi*wh/max_tau;
    start_a = mean(y(round(max_tau/2):end));
    start_c = y(1);
    
    fo = fitoptions('Method','NonLinearLeastSquares','Lower',[0, 0, 0],'Upper',[Inf, Inf, 2*y(1)],'StartPoint',[start_a, start_b, start_c]);%,'Weights',1./err_y.^2);   % options for the fit
    fit_out = fit(tau,y,ft,fo);
    A(q_c) = fit_out.a;             % saturation value
    qvel(q_c) = fit_out.b;      % q*v (this is what this fit is done for)
    B(q_c) = fit_out.c;             % offset
    temp = 0.5*diff(confint(fit_out,0.682));    % this and next three lines: errors on A, tau_q and B
    err_A(q_c) = temp(1);
    err_qvel(q_c) = temp(2);
    err_B(q_c) = temp(3);
    
    %     fit_out = fminsearch('funz_fitbal', [start_a, start_b, start_c],optimset('MaxFunEvals',1e7,'MaxIter',1e7,'TolFun',1e-10), y,  tau(1:end-q_c+1), 1./err_y.^2);
    %     A(q_c) = fit_out(1);
    %     qvel(q_c) = fit_out(2);
    %     B(q_c) = fit_out(3);
    
%     chi2_ballistic(q_c) = sum( (y - A(q_c)*(1-sin(qvel(q_c)*tau)./(qvel(q_c)*tau))+B(q_c)).^2 ./ err_y.^2 ) / (max_tau-4);
    
    % plot stuff
    cla;
    semilogx(tau,y,'.');
    %     errorbar(tau,y,err_y,'.');
    hold on;
    semilogx(tau, A(q_c)*(1-sin(qvel(q_c)*tau)./(qvel(q_c)*tau))+B(q_c), 'r');
    title(['q = ',num2str(q_c)],'FontSize',14);
    xlabel('Lag time \tau [s]', 'FontSize',14);
    ylabel('Intensity [a.u.]', 'FontSize',14);
    shg
    pause(0.01);
    
end


%% fit of (q*v) vs q

% change the next two lines if you want to fit a different range (in q) of the datapoints
min_q_fit = 1;
max_q_fit = length(q);

% plot stuff
clf;
hold on;
% errorbar(q(min_q_fit:max_q_fit),tau_q(min_q_fit:max_q_fit),err_tau_q(min_q_fit:max_q_fit),'g.');
dummy_q = [min_q_fit:max_q_fit]';
errorbar(dummy_q(min_q_fit:max_q_fit),qvel(min_q_fit:max_q_fit),err_qvel(min_q_fit:max_q_fit),'g.');
ylim([ min(qvel)-1 max(qvel)+1]);
%set(gca,'XScale','log','YScale','log');

% interactive stuff
disp('Click on the leftmost and on the rightmost of the points you want to fit, then press enter');
title('Click on the leftmost and on the rightmost of the points you want to fit, then press enter');
[x,y] = getpts(gcf);
if x(2)<x(1)
    x= x (end:-1:1);
    y = y(end:-1:1);
end
% [~,min_q_fit]=min((x(1)-q).^2+(y(1)-tau_q).^2);
% [~,max_q_fit]=min((x(2)-q).^2+(y(2)-tau_q).^2);
[~,min_q_fit]=min((x(1)-dummy_q).^2+(y(1)-qvel).^2);
[~,max_q_fit]=min((x(2)-dummy_q).^2+(y(2)-qvel).^2);


%%
% plot stuff
clf;
hold on;
errorbar(q,qvel,err_qvel,'b.');
ylim([ min(qvel)-1 max(qvel)+1]);

qq = q(min_q_fit:max_q_fit);
qqvel = qvel(min_q_fit:max_q_fit);
err_qqvel = err_qvel(min_q_fit:max_q_fit);
errorbar(qq,qqvel,err_qqvel,'g.');

disp('Drag a rectangle around outliers, then press enter');
title('Drag a rectangle around outliers,, then press enter');
outliers_box = getrect(gca);

outliers = qq > outliers_box(1) & qq < outliers_box(1)+outliers_box(3) & qqvel > outliers_box(2) & qqvel < outliers_box(2) + outliers_box(4);
errorbar(qq(outliers),qqvel(outliers),err_qqvel(outliers),'r.');

qq(outliers) = [];
qqvel(outliers) = [];
err_qqvel(outliers) = [];


%set(gca,'XScale','log','YScale','log');
xlabel('q [px^{-1}]','FontSize',14);
ylabel('\tau_q [s]','FontSize',14);


% actual fit



foo2 = fit(qq,qqvel,fittype('vel*x+c','independent','x'),fitoptions('Method','NonLinearLeastSquares','Weights',1./err_qqvel.^2));

% plot stuff
plot(qq,foo2.vel*qq+foo2.c,'k');
% set(gca,'XScale','log','YScale','log');
xlabel('q [px^{-1}]','FontSize',14);
ylabel('\tau_q [s]','FontSize',14);
% conversion to mum2/s
estimated_vel = foo2.vel * px2mum;       % now in mum2/s (foo.Dm is in px2/s)
disp(['velocity = ', num2str(estimated_vel)]);

% save(filename);