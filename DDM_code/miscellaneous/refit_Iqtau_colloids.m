
scroll_through_Iqtau_fits(1e4.*Iqtau,fits_array,tau);

%% fit of the time constant tau_q vs q

% change the next two lines if you want to fit a different range (in q) of the datapoints
min_q_fit = 1;
max_q_fit = length(q);

% plot stuff
figure(103);
clf;
hold on;
% errorbar(q(min_q_fit:max_q_fit),tau_q(min_q_fit:max_q_fit),err_tau_q(min_q_fit:max_q_fit),'g.');
% dummy_q = [min_q_fit:max_q_fit]';                                                     %use this line instead of the following if you get an unexpected error
dummy_q = [1:length(q)]';
errorbar(dummy_q(min_q_fit:max_q_fit),tau_q(min_q_fit:max_q_fit),err_tau_q(min_q_fit:max_q_fit),'g.');
plot(dummy_q(min_q_fit:max_q_fit),tau_q(min_q_fit:max_q_fit),'+','MarkerSize',1);
set(gca,'XScale','log','YScale','log');
% xlabel('q [px^{-1}]','FontSize',14);
xlabel('mode','FontSize',14);
ylabel('\tau_q [s]','FontSize',14);

% interactive stuff
disp('Click on the leftmost and on the rightmost of the points you want to fit, then press enter');
title('Click on the leftmost and on the rightmost of the points you want to fit, then press enter');
[x,y] = getpts(gcf);
if x(2)<x(1)
    x = x(end:-1:1);
    y = y(end:-1:1);
end

% disp(['(', num2str(x(1)), ', ', num2str(y(1)), ')']);
% disp(['(', num2str(x(2)), ', ', num2str(y(2)), ')']);

% [~,min_q_fit]=min((x(1)-q).^2+(y(1)-tau_q).^2);
% [~,max_q_fit]=min((x(2)-q).^2+(y(2)-tau_q).^2);
[~,min_q_fit]=min((x(1)-dummy_q).^2+(y(1)-tau_q).^2);
[~,max_q_fit]=min((x(2)-dummy_q).^2+(y(2)-tau_q).^2);

% disp(['Min_q_fit = ', num2str(min_q_fit)]);
% disp(['Max_q_fit = ', num2str(max_q_fit)]);

%%
% plot stuff
figure(103);
clf;
hold on;
errorbar(q([1:min_q_fit,max_q_fit:end]),tau_q([1:min_q_fit,max_q_fit:end]),err_tau_q([1:min_q_fit,max_q_fit:end]),'b.');
errorbar(q(min_q_fit:max_q_fit),tau_q(min_q_fit:max_q_fit),err_tau_q(min_q_fit:max_q_fit),'g.');

set(gca,'XScale','log','YScale','log');
xlabel('q [px^{-1}]','FontSize',14);
ylabel('\tau_q [s]','FontSize',14);


% actual fit
foo = fit(q(min_q_fit:max_q_fit),tau_q(min_q_fit:max_q_fit),fittype('1/(Dm*x^2)','independent','x'),fitoptions('Method','NonLinearLeastSquares','Weights',1./err_tau_q(min_q_fit:max_q_fit).^2));
foo2 = fit(q(min_q_fit:max_q_fit),tau_q(min_q_fit:max_q_fit),fittype('1/(Dm*x^alpha)','independent','x'),fitoptions('Method','NonLinearLeastSquares','StartPoint',[1, 2],'Weights',1./err_tau_q(min_q_fit:max_q_fit).^2));

% plot stuff
plot(q(min_q_fit:max_q_fit),1./(foo.Dm*q(min_q_fit:max_q_fit).^2),'r');
plot(q(min_q_fit:max_q_fit),1./(foo2.Dm*q(min_q_fit:max_q_fit).^foo2.alpha),'k');

set(gca,'XScale','log','YScale','log');
xlabel('q [px^{-1}]','FontSize',14);
ylabel('\tau_q [s]','FontSize',14);
legend('non-fitted data points','fitted data points','fit with exponent fixed at 2',['fit with free exponent \alpha = ',num2str(foo2.alpha)],'Location','SouthWest')
% conversion to mum2/s
Dm = foo.Dm * px2mum^2;       % now in mum2/s (foo.Dm is in px2/s)
DmFo = foo2.Dm * px2mum^2;
disp(['Diffusion Coefficient (exponent fixed at 2) = ', num2str(Dm),' um2/s']);
disp(['Diffusion Coefficient (free exponent) = ', num2str(DmFo),' um2/s']);

