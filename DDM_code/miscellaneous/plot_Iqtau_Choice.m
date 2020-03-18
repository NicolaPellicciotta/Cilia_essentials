function [ min_q_fit, max_q_fit ] = plot_Iqtau_Choice(q, tau_q_array, index)
%% Plots the data passed to it and extracts revised fitting range interactively

%%
min_q_fit = 1;     % Define range in q of datapoints
max_q_fit = length(q);
%a = 'min_q_fit:max_q_fit'; %b = '1:min_q_fit'; %c = 'max_q_fit:end';
dummy_q = (min_q_fit:max_q_fit)';      %use this line instead of the following if you get an unexpected error dummy_q = [1:length(q)]';

%% Plot figure
% Which fit?
switch index
    case 1
        yla = '\tau_q / s';
        Sca = 'log';
    case 3
        yla = '\omega /s';
        Sca = 'log';
end

% Axes
clf;        
hold on;
set(gca,'XScale',Sca,'YScale',Sca);
xlabel('mode','FontSize',14);           % c.f. % xlabel('q [px^{-1}]','FontSize',14);
ylabel(yla,'FontSize',14);

% Plots
errorbar(dummy_q(min_q_fit:max_q_fit),tau_q_array(min_q_fit:max_q_fit,index), ...   % index 1 2 3 4
            tau_q_array(min_q_fit:max_q_fit,(index+1)),'g.');                       % tau_q_array b berr (d derr)
plot(dummy_q(min_q_fit:max_q_fit),tau_q_array(min_q_fit:max_q_fit,index),'+','MarkerSize',1);

% Interactive part
Request = ['Click on the leftmost and on the rightmost ' ...
            'of the points you want to fit, then press enter'];
disp(Request);
title(Request);
[x,y] = getpts(gcf);    
if x(2)<x(1)            % Does not matter if 'wrong' way round
    x = x(end:-1:1);
    y = y(end:-1:1);
end                     

[~,min_q_fit]=min((x(1)-dummy_q).^2+(y(1)-tau_q_array(:,1)).^2); % Find nearest q modes ...
[~,max_q_fit]=min((x(2)-dummy_q).^2+(y(2)-tau_q_array(:,1)).^2);   
disp(['Min_q_fit = ', num2str(min_q_fit)]);                      % ... and show them
disp(['Max_q_fit = ', num2str(max_q_fit)]);

end