%% this script is meant to be used after MultiDDM analysis on many samples 
%% with the AccumData.mat, see MultiDDM routines for more info 


clearvars -except dc
%   load the AccumData variable 
load('/home/np451/Desktop/Odiri/run1/AccumData.mat',...
    'SampleType',...
    'figures_folder',..._nosparse
    'sampletypes',...
    'timepoints',...
    'donors',...
    'inserts',...
    'positions',...
    'boxsizes_vector',...
    'q_limits_1oum');


for i = 1:numel(sampletypes)
    for j = 1:numel(timepoints)
        
        MergedData_alldonors(i, j) = merge_SampleType_data(SampleType,...
                i, j, donors, inserts, '');
        
        % separate donors: we have 3 donors, 4 drugs, 4 timepoints
        for dc = 1:numel(donors)
            
            MergedData(dc,i,j) = merge_SampleType_data(SampleType,...
                i, j, donors{dc}, inserts, '');
            
        end %for dc
        
    end %for j
end %for i

figures_folder = '/home/np451/Desktop/Odiri/run1/figures/';
%% now figure of histograms (from this pick a fig 4A, B, C, D)

darkB = [0 0.4470 0.7410];
darkO = [0.8500 0.3250 0.0980];
darkG = [0.13 0.53 0.13];

lightB = darkB + 0.25;
lightO = darkO + [0.15 0.23 0.2];
darkGrey = [0.15 0.15 0.15];
lightGrey = [0.65 0.65 0.65];
cmap = [27,158,119 ;
    217,95,2 ;
    117,112,179 ] ./ 255;

st_cmap = vertcat(darkB, darkG, darkO, darkGrey) ;

%%close all


for i = 1:numel(MergedData) 
    
    [~,stc,~] = ind2sub(  size(MergedData),i);
    
    
    lfs = 10;   % labels fontsize
    tfs = 8;    % ticks fontsize
    binedges = 0:0.5:30;
%     binedges = 0:0.4:15;
    
    bsi = find(boxsizes_vector == 64);
    
    hf = figure;
    hf.Units = 'centimeters';
    hf.Position([3 4]) = [8 6];
    
    ha = axes;
    ha.Box = 'on';
    ha.NextPlot = 'add';
    ha.FontSize = tfs;
    ha.FontName = 'Arial';
    
    
%     histogram(MergedData(i).Frequency_Hz{bsi}, binedges,...
%         'Normalization','Probability',...
%         'DisplayStyle','Stairs',...
%         'EdgeColor', darkB,...
%         'LineStyle', '-',...
%         'LineWidth',1.5);
    histogram(MergedData(i).Frequency_Hz{bsi}, binedges,...
        'Normalization','Probability',...
        'FaceColor', st_cmap(stc,:),...
        'FaceAlpha', 1,...
        'EdgeColor', 'none');
    
    ha.XLabel.String = 'CBF, [Hz]';
    ha.XLabel.FontSize = lfs;
    
    ha.YLabel.String = 'Counts, normalised';
    ha.YLabel.FontSize = lfs;
    
    % fix title
    ha.Title.String = strjoin({MergedData(i).donors_str{:},MergedData(i).timepoint_str},' ');
    ha.Title.FontWeight = 'normal';
    ha.Title.FontSize = lfs;
    
    ha.YLim = [0 0.45];
    ha.XLim = [0 20];

    hleg = legend({MergedData(i).sampletype_str});
    hleg.String = replace(hleg.String, {'LumacaFTOR';'Orkambi';'T1a'},{'Lumacaftor';'Lumacaftor + Ivacaftor';'T\alpha1'});

    hleg.Box = 'off';
    hleg.FontSize = lfs;
%     hleg.Interpreter = 'none';

    hleg.Position(1) = sum(ha.Position([1 3])) - hleg.Position(3);
    hleg.Position(2) = sum(ha.Position([2 4])) - 1.3*hleg.Position(4);
    
    savename = strjoin({'CBF_histogram',...
        MergedData(i).donors_str{:},...
        MergedData(i).sampletype_str,...
        MergedData(i).timepoint_str},'_');
    savename = fullfile(figures_folder,savename);
   % print2svg(hf, savename, 1, 1)
    
end

%% plot each donors sigmoids on the same set of axes (one per donor) 
% so 5 sigmoids per set of axes


%%close all

markarr = '^odsv';
linest = {'-','-','--','--','-.'};


darkB = [0 0.4470 0.7410];
darkO = [0.8500 0.3250 0.0980];
darkG = [0.13 0.53 0.13];

lightB = darkB + 0.25;
lightO = darkO + [0.15 0.23 0.2];
darkGrey = [0.15 0.15 0.15];
lightGrey = [0.65 0.65 0.65];
cmap = [27,158,119 ;
    217,95,2 ;
    117,112,179 ] ./ 255;

st_cmap = vertcat(darkB, darkG, darkO, darkGrey) ;


lfs = 10;   % labels fontsize
tfs = 8;    % ticks fontsize


for stc = 1:numel(sampletypes)
    for tpc = 1:numel(timepoints)
        
        
        hf = figure;
        ha = axes;
        ha.YLim = [0 8];
        [hf, ha, he, hpf, hpv, hpp] = plot_sigmoid_errorbar_compare(MergedData(:,stc, tpc), 'lines', ha, 1,1,'left');
        
        
        hf.Units = 'centimeters';
        hf.Position([3 4]) = [8 6];
        
        ha.FontSize = tfs;
        ha.FontName = 'Arial';
        
        ha.XLabel.FontSize = lfs;
        ha.YLabel.FontSize = lfs;
        
        ha.YLim = [0 9];
        
        for dc = 1:numel(donors)
            
            he(dc).Color = cmap(dc,:);
            he(dc).MarkerEdgeColor = cmap(dc,:);
            hpf(dc).Color = cmap(dc,:);
            hpv(dc).Color = cmap(dc,:);
            hpp(dc).FaceColor = cmap(dc,:);
            
            he(dc).Marker = markarr(dc);
            %         hpv(i).LineStyle = linest{i};
            %         hpf(i).LineStyle = linest{i};
        end %for dc
        
        uistack(hpv,'bottom');
        uistack(hpp,'bottom');
        %     uistack(he,'top');
        
        hl = findobj(hf.Children,'type','legend');
        hl.FontSize = lfs;
        hl.Units = 'normalized';
        hl.Position(1) = ha.Position(1);
        hl.Position(2) = sum(ha.Position([2 4])) - hl.Position(4) - 0.01;
        
        % fix title
        ha.Title.String = strjoin({sampletypes{stc}, timepoints{tpc} },' ');
        ha.Title.String = replace(ha.Title.String, {'LumacaFTOR';'Orkambi';'T1a'},{'Lumacaftor';'Lumacaftor + Ivacaftor';'T\alpha1'});
        ha.Title.FontWeight = 'normal';
        ha.Title.FontSize = lfs;
        
        
        drawnow
        
        savename = strjoin({'Sigmoids',...
            sampletypes{stc},...
            timepoints{tpc}},'_');
        disp(savename)
        savename = fullfile(figures_folder,savename);
       % print2svg(hf, savename, 1, 1);
    end %for tpc
end %for stc

%% plot sigmoids shoulders vs time, first for each donor, then as geometric mean/std of the donors' shoulders. 


% %%close all

markarr = '^odsv';
linest = {'-','-','--','--','-.'};


darkB = [0 0.4470 0.7410];
darkO = [0.8500 0.3250 0.0980];
lightB = darkB + 0.25;
lightO = darkO + [0.15 0.23 0.2];
darkG = [0.13 0.53 0.13];
darkGrey = [0.15 0.15 0.15];
lightGrey = [0.65 0.65 0.65];
cmap = [27,158,119 ;
    217,95,2 ;
    117,112,179 ] ./ 255;


lfs = 10;   % labels fontsize
tfs = 8;    % ticks fontsize


flag_normalised = 0;
whichshoulder = 'left';


% if flag_normalised
%     dampingname = 'Damping_1ocycles_fit_out';
%     parname = 'b';
% else
%     dampingname = 'Damping_Hz_fit_out2';
%     parname = 'mu';
% end

switch whichshoulder
    case 'left'
        lambda2_mult = exp(-2);
        lambda2_shift = -2/log(10);
    case 'right'
        lambda2_mult = exp(2);
        lambda2_shift = 2/log(10);
    otherwise
        error('either ''left'' or ''right''.');
end


% save where the shoulder is for all sigmoids
shiftedmuMat = nan(size(MergedData));
confint_shiftedmuMat = nan(size(MergedData)); % half the confidence interval


for i = 1:numel(MergedData)

    % read inflexion point with shift (this is in Log scale)
    shiftedmuMat(i) = MergedData(i).Damping_Hz_fit_out2.mu + lambda2_shift;
    
    confint_shiftedmuMat(i) = diff(par_confint(MergedData(i).Damping_Hz_fit_out2,'mu',0.68))./2 ;
    
    timeMat(i) = str2double(MergedData(i).timepoint_str(1:2));
        
end

timeMat(isnan(timeMat))=0; %%%%% correction for 0h
shiftedmuMat = reshape(shiftedmuMat,size(MergedData));
confint_shiftedmuMat = reshape(confint_shiftedmuMat,size(MergedData));
timeMat = reshape(timeMat,size(MergedData));

%this gives me error
%plottimeMat = 1+timeMat./24 + repmat(linspace(-0.1, 0.1, 4),3,1,4);
plottimeMat = 1+timeMat./48 + repmat(linspace(-0.1, 0.1, 3),3,1,2);

st_cmap = vertcat(darkB, darkG, darkO, darkGrey) ;

flag_no72h = 0;

%%%%% here I am trying to understand I think it is the timepoints
if flag_no72h
    ind_plot = 1:3;
else
    ind_plot = 1:numel(timepoints);
end

% figure with all donors on same axes
lfs = 10;   % labels fontsize
tfs = 8;    % ticks fontsize

hf = figure;
hf.Units = 'centimeters';
hf.Position([3 4]) = [8 6];

ha = axes;
ha.NextPlot = 'add';
ha.Box = 'on';
ha.FontName = 'Arial';
ha.FontSize = tfs;

clear he;

% stagger time for each donor
for dc = 1:numel(donors)
    
    
    % plot data
    he(dc,:) = plot(squeeze(plottimeMat(dc,:,ind_plot))',...
        10.^squeeze(shiftedmuMat(dc,:,ind_plot))');
      
    
    % appearances
    for stc = 1:numel(sampletypes)
    he(dc,stc).Marker = markarr(dc);
    he(dc,stc).MarkerFaceColor = 'w';
    he(dc,stc).MarkerSize = 6;
    %     he_w(dc).Color = lightO;
    he(dc,stc).Color = st_cmap(stc,:);
    he(dc,stc).LineWidth = 1.2;
%     he(dc).CapSize = 5;
    he(dc,stc).LineStyle = 'none';
    end
    
end

legend
hl = legend(he(1:3:end),{MergedData(1,:,1).sampletype_str},{MergedData(1,:,1).sampletype_str})
 %   replace({MergedData(1,:,1).sampletype_str}, {'LumacaFTOR';'Orkambi';'T1a'},{'Lumacaftor';'Lumacaftor + Ivacaftor';'T\alpha1'}));


hl.Box = 'on';
hl.Orientation = 'vertical';
hl.FontSize = lfs;
hl.Units = 'normalized';

% axes properties
ha.XTick = 1:size(MergedData,2);
%ha.XLim = minmax(ind_plot) + 0.5.*[-1 1];
ha.XLim = [0.5, 2.8];

ha.XTickLabel = {MergedData(1,1,ind_plot).timepoint_str} ;
% ha.XAxis.FontSize = lfs;
ha.XLabel.String = 'Time';
ha.XLabel.FontSize = lfs;
ha.YLabel.String = '\lambda^2, [\mum^2]';
ha.YLabel.FontSize = lfs;
ha.YLim = [1e2 10^3.1];
setsemilogy
% 
% positions
ha.Position = [0.1591 0.1546 0.82 0.79];
% hl.Position(1) = ha.Position(1) - 0.02;
% hl.Position(2) = sum(ha.Position([2 4])) - hl.Position(4);
ha.YLabel.Units = 'normalized';
ha.YLabel.Position(1) = -0.13;
ha.YLabel.VerticalAlignment = 'baseline';

savename = fullfile('left_lambda2(t)');
if flag_no72h
    savename = [savename,'_no72h'];
end
savename = fullfile(figures_folder,savename);
%print2svg(hf, savename, 1, 1);

%% Now only plot of averaged quantities

% weighted mean
[shiftedmu_donavg, shiftedmu_err] = weightednanmean(shiftedmuMat,confint_shiftedmuMat,1);
shiftedmu_donavg = squeeze(shiftedmu_donavg);
shiftedmu_err = squeeze(shiftedmu_err);
Lambda2Mat_donavg = 10.^shiftedmu_donavg;
Lambda2Mat_uperr = 10.^(shiftedmu_donavg + shiftedmu_err) - Lambda2Mat_donavg ;
Lambda2Mat_loerr =Lambda2Mat_donavg - 10.^(shiftedmu_donavg - shiftedmu_err);

timeMat_donavg = squeeze(mean(plottimeMat,1));
% now first dimension is treatment, second is timepoint

%%close all
flag_trendline = 1;
flag_no72h = 0;

if flag_no72h
    ind_plot = 1:3;
else
    ind_plot = 1:2;
end


lfs = 10;   % labels fontsize
tfs = 8;    % ticks fontsize

% plot errorbar

hf = figure;
hf.Units = 'centimeters';
hf.Position([3 4]) = [8 6];

ha = axes;
ha.NextPlot = 'add';
ha.Box = 'on';
ha.FontName = 'Arial';

% plot data
he = errorbar(timeMat_donavg(:,ind_plot)',...
    Lambda2Mat_donavg(:,ind_plot)',...
    Lambda2Mat_loerr(:,ind_plot)',...
    Lambda2Mat_uperr(:,ind_plot)',...
    'o');


% legend
hl = legend(he, replace({MergedData(1,:,1).sampletype_str},...
    {'LumacaFTOR';'Orkambi';'T1a'},{'Lumacaftor';'Lumacaftor + Ivacaftor';'T\alpha1'}));

for stc = 1:numel(sampletypes)
    
    % appearances
    he(stc).Marker = 'o';
    he(stc).MarkerFaceColor = 'w';
    he(stc).LineWidth = 1.2;
    he(stc).CapSize = 5;
    % set(he_m,'LineStyle','none');
    he(stc).Color = st_cmap(stc,:);
    
    % create trendline
    if flag_trendline
        fitline = fit(mean(timeMat_donavg(:,ind_plot))', shiftedmu_donavg(stc,ind_plot)','poly1',...
            fitoptions('weights',1./shiftedmu_err(stc,ind_plot)));
        
        hpf(stc) = plot(timeMat_donavg(stc,ind_plot), 10.^fitline(timeMat_donavg(stc,ind_plot)));
        hpf(stc).Marker = 'none';
        hpf(stc).LineWidth = 1.2;
        hpf(stc).Color = st_cmap(stc,:);
    end
    
end

hl.Box = 'off';
hl.FontSize = lfs;
hl.Units = 'normalized';

% axes properties
ha.FontSize = 8;
ha.XTick = ind_plot;
ha.XLim = minmax(ha.XTick) + 0.5.*[-1 1];
ha.XTickLabel = {MergedData(1,1,ind_plot).timepoint_str} ;
% ha.XAxis.FontSize = lfs;
ha.XLabel.String = 'Time';
ha.XLabel.FontSize = lfs;
ha.YLabel.String = '<\lambda^2>_d, [\mum^2]';
ha.YLabel.FontSize = lfs;
ha.YLim = [1e2 1e5];
setsemilogy

% positions
ha.Position = [0.1591 0.1546 0.82 0.79];
% hl.Position(1) = ha.Position(1) - 0.01;
% hl.Position(2) = ha.Position(2) + 0.01;
hl.Position(1) = ha.Position(1) - 0.02;
hl.Position(2) = sum(ha.Position([2 4])) - hl.Position(4);
ha.YLabel.Units = 'normalized';
ha.YLabel.Position(1) = -0.13;
ha.YLabel.VerticalAlignment = 'baseline';

savename = 'left_lambda2(t)_average';
if flag_trendline
    savename = [savename,'_trendline'];
end

if flag_no72h
    savename = [savename,'_no72h'];
end
savename = fullfile(figures_folder,savename);
%print2svg(hf, savename, 1, 1);

%% now plot log10(lambda2(48h)) - log10(lambda2(0h)) - [ log10(lambda2(48h)) - log10(lambda2(0h))]ctrl
% for x donor treatment but with donor on the x axis
%%close all

% %%close all

markarr = '^odsv<>';
linest = {'-','-','--','--','-.'};


darkB = [0 0.4470 0.7410];
darkO = [0.8500 0.3250 0.0980];
lightB = darkB + 0.25;
lightO = darkO + [0.15 0.23 0.2];
darkG = [0.13 0.53 0.13];
darkGrey = [0.15 0.15 0.15];
lightGrey = [0.65 0.65 0.65];


cmap = [166,206,227
    31,120,180
    178,223,138] ./ 255;

st_cmap = [27,158,119
    217,95,2
    152,67,1
    117,112,179
    231,41,138
    102,166,30 ] ./ 255;
st_cmap = vertcat(st_cmap, darkGrey);

st_cmap=[1 0 0
    1 0 0
    0 0 1];
x_stagger = [-0.09 0 0.09];

lfs = 10;   % labels fontsize
tfs = 8;    % ticks fontsize
lw = 1.2;

%[hf, ha] = figure_template;
hf = figure;
ha = axes;

clear he;

matforpvalues = nan(3,4);
errmatforpvalues = nan(3,4);
delta2_lambda2_48h_mat = []; %y to plot
le_delta2_lambda2_48h_mat = []; %yneg
ue_delta2_lambda2_48h_mat = []; %ypos

% sampletype on x axes, 
for stc = 2:numel(sampletypes)
    for dc = 1:numel(donors)
        
        % no need to apply the left shift as we are then subtracting one
        % from the other
        mu_48h = MergedData(dc,stc,2).Damping_Hz_fit_out2.mu;
        mu_48h_ci = diff(par_confint(MergedData(dc,stc,2).Damping_Hz_fit_out2,'mu',0.68))./2;
        mu_00h = MergedData(dc,stc,1).Damping_Hz_fit_out2.mu;
        mu_00h_ci = diff(par_confint(MergedData(dc,stc,1).Damping_Hz_fit_out2,'mu',0.68))./2;
        
        delta_mu = mu_48h - mu_00h;
        err_delta_mu = sqrt(mu_48h_ci.^2 + mu_00h_ci.^2);
        
        % now same for control (yeah I shouldn;t recalculate it every time
        % but whatever)
        mu_48h_ctrl = MergedData(dc,1,2).Damping_Hz_fit_out2.mu;
        mu_48h_ci_ctrl = diff(par_confint(MergedData(dc,1,2).Damping_Hz_fit_out2,'mu',0.68))./2;
        mu_00h_ctrl = MergedData(dc,1,1).Damping_Hz_fit_out2.mu;
        mu_00h_ci_ctrl = diff(par_confint(MergedData(dc,1,1).Damping_Hz_fit_out2,'mu',0.68))./2;
        
        delta_mu_ctrl = mu_48h_ctrl - mu_00h_ctrl;
        err_delta_mu_ctrl = sqrt(mu_48h_ci_ctrl.^2 + mu_00h_ci_ctrl.^2);
        
        delta2_mu = delta_mu - delta_mu_ctrl;
        err_delta2_mu = sqrt(err_delta_mu.^2 + err_delta_mu_ctrl.^2);
        
        delta2_lambda2 = 10.^delta2_mu;
        10.^delta_mu ./ 10.^delta_mu_ctrl;
        ue_delta2_lambda2 = 10.^(delta2_mu + err_delta2_mu) - 10.^(delta2_mu);
        le_delta2_lambda2 = 10.^(delta2_mu) - 10.^(delta2_mu - err_delta2_mu);
                
        delta2_lambda2_48h_mat(dc,stc) = delta2_lambda2; %y to plot
        le_delta2_lambda2_48h_mat(dc,stc) = le_delta2_lambda2; %yneg
        ue_delta2_lambda2_48h_mat(dc,stc) = ue_delta2_lambda2; %ypos
        
%         he(stc,dc) = errorbar(stc + x_stagger(dc), delta_mu, err_delta_mu);
 %       he(stc,dc) = errorbar(stc + x_stagger(dc), delta2_lambda2,...
       he(stc,dc) = errorbar(dc + x_stagger(stc), delta2_lambda2,...
         le_delta2_lambda2, ue_delta2_lambda2);
hold on;
        
        he(stc,dc).Marker = markarr(dc);
        he(stc,dc).MarkerFaceColor = 'w';
        he(stc,dc).Color = st_cmap(stc,:);
 %       he(stc,dc).Color = st_cmap(stc);

         he(stc,dc).CapSize = 6;
         
        he(stc,dc).LineWidth = lw;
        he(stc,dc).LineStyle = 'none';
 %       he(stc,dc).DisplayName = MergedData(dc,stc,1).donors_str{:};
               he(stc,dc).DisplayName = MergedData(dc,stc,1).sampletype_str;
        
        matforpvalues(dc,stc) = delta2_mu;
        errmatforpvalues(dc,stc) = 2*err_delta2_mu; %twice because err_delta_mu was halved for plotting
        
    end %for dc
end %for

save('lambda2_ratio(48h)_CtrlNormalised.mat',...
    'delta2_lambda2_48h_mat',...
    'le_delta2_lambda2_48h_mat',...
    'ue_delta2_lambda2_48h_mat',...
    'sampletypes');

setsemilogy




% axes properties
%ha.XTick = 1:numel(sampletypes)-1;
ha.XTick = 1:numel(donors);

%ha.XLim = minmax(ha.XTick) + 0.5.*[-1 1];
ha.XLim = [0.8,3.2];%minmax(ha.XTick) + 0.5.*[-1 1];
ha.YLim = [1e-2,1e2];
ha.XTickLabel = {'CF1','CF2','CF3'};
ha.XTickLabelRotation = 15;

plot(xlim,ones(2,1),'--','color',[0.5 0.5 0.5],'linewidth',0.4)

% ha.XLabel.String = ' ';
% ha.XLabel.FontSize = lfs;
ha.YLabel.String = '(\lambda^2_{48h} / \lambda^2_{0h}) / (\lambda^2_{48h} / \lambda^2_{0h})_{ctrl}';
% ht = fix_xticklabels(ha);
% set(ht,'FontSize',lfs,...
%     'FontName','Arial',...
%     'Color',darkGrey);

hl = legend(he(2:end,end));
hl.String = replace(hl.String, {'_L';'_R'},{'Lumacaftor + Ivacaftor';'RPL554'});
hl.Box = 'on';
hl.FontSize = tfs*2;
hl.Units = 'normalized';
hl.String = regexprep(hl.String,'d(?=\d)','s');
% positions
%ha.Position = hapos + [0.03 0 0 0];
%hl.Position(1) = ha.Position(1) - 0.05;
%hl.Position(2) = ha.Position(2) + ha.Position(4) - hl.Position(4) - 0.01;
hl.Location = 'northeast'
ha.YLabel.Units = 'normalized';
ha.YLabel.Position(1) = -0.13;
ha.YLabel.VerticalAlignment = 'baseline';


savename = 'lambda2_ratio_48h_normalised_by_control_donors.pdf';
%print2svg(hf, savename, 1);
%print2pdf(hf, savename, 0);
saveas(gcf,savename)
% 
% disp('using paired t-test')
% for stc = 2:numel(sampletypes)
%     disp(sampletypes{stc})
%     disp('P value \lambda2(48h) / \lambda2(0h) vs control');
%     [~,p] = ttest(matforpvalues(:,stc), 0,'tail','left');
%     disp(p)
% end



matforpvalues
errmatforpvalues

%%

%%%% now plot log10(lambda2(48h)) - log10(lambda2(0h)) - [ log10(lambda2(48h)) - log10(lambda2(0h))]ctrl
% for each donor/treatment

%%close all

% %%close all

markarr = '^odsv<>';
linest = {'-','-','--','--','-.'};


darkB = [0 0.4470 0.7410];
darkO = [0.8500 0.3250 0.0980];
lightB = darkB + 0.25;
lightO = darkO + [0.15 0.23 0.2];
darkG = [0.13 0.53 0.13];
darkGrey = [0.15 0.15 0.15];
lightGrey = [0.65 0.65 0.65];


cmap = [166,206,227
    31,120,180
    178,223,138] ./ 255;

st_cmap = [27,158,119
    217,95,2
    152,67,1
    117,112,179
    231,41,138
    102,166,30 ] ./ 255;
st_cmap = vertcat(st_cmap, darkGrey);

st_cmap=[1 0 0
    1 0 0
    0 0 1];
x_stagger = [-0.09 0 0.09];

lfs = 10;   % labels fontsize
tfs = 8;    % ticks fontsize
lw = 1.2;

%[hf, ha] = figure_template;
hf = figure;
ha = axes;

clear he;

matforpvalues = nan(3,4);
errmatforpvalues = nan(3,4);
delta2_lambda2_48h_mat = []; %y to plot
le_delta2_lambda2_48h_mat = []; %yneg
ue_delta2_lambda2_48h_mat = []; %ypos

% sampletype on x axes, 
for stc = 2:numel(sampletypes)
    for dc = 1:numel(donors)
        
        % no need to apply the left shift as we are then subtracting one
        % from the other
        mu_48h = MergedData(dc,stc,2).Damping_Hz_fit_out2.mu;
        mu_48h_ci = diff(par_confint(MergedData(dc,stc,2).Damping_Hz_fit_out2,'mu',0.68))./2;
        mu_00h = MergedData(dc,stc,1).Damping_Hz_fit_out2.mu;
        mu_00h_ci = diff(par_confint(MergedData(dc,stc,1).Damping_Hz_fit_out2,'mu',0.68))./2;
        
        delta_mu = mu_48h - mu_00h;
        err_delta_mu = sqrt(mu_48h_ci.^2 + mu_00h_ci.^2);
        
        % now same for control (yeah I shouldn;t recalculate it every time
        % but whatever)
        mu_48h_ctrl = MergedData(dc,1,2).Damping_Hz_fit_out2.mu;
        mu_48h_ci_ctrl = diff(par_confint(MergedData(dc,1,2).Damping_Hz_fit_out2,'mu',0.68))./2;
        mu_00h_ctrl = MergedData(dc,1,1).Damping_Hz_fit_out2.mu;
        mu_00h_ci_ctrl = diff(par_confint(MergedData(dc,1,1).Damping_Hz_fit_out2,'mu',0.68))./2;
        
        delta_mu_ctrl = mu_48h_ctrl - mu_00h_ctrl;
        err_delta_mu_ctrl = sqrt(mu_48h_ci_ctrl.^2 + mu_00h_ci_ctrl.^2);
        
        delta2_mu = delta_mu - delta_mu_ctrl;
        err_delta2_mu = sqrt(err_delta_mu.^2 + err_delta_mu_ctrl.^2);
        
        delta2_lambda2 = 10.^delta2_mu;
        10.^delta_mu ./ 10.^delta_mu_ctrl;
        ue_delta2_lambda2 = 10.^(delta2_mu + err_delta2_mu) - 10.^(delta2_mu);
        le_delta2_lambda2 = 10.^(delta2_mu) - 10.^(delta2_mu - err_delta2_mu);
                
        delta2_lambda2_48h_mat(dc,stc) = delta2_lambda2; %y to plot
        le_delta2_lambda2_48h_mat(dc,stc) = le_delta2_lambda2; %yneg
        ue_delta2_lambda2_48h_mat(dc,stc) = ue_delta2_lambda2; %ypos
        
%         he(stc,dc) = errorbar(stc + x_stagger(dc), delta_mu, err_delta_mu);
        he(stc,dc) = errorbar(stc + x_stagger(dc), delta2_lambda2,...
            le_delta2_lambda2, ue_delta2_lambda2);
hold on;
        
        he(stc,dc).Marker = markarr(dc);
        he(stc,dc).MarkerFaceColor = 'w';
        he(stc,dc).Color = st_cmap(stc,:);
 %       he(stc,dc).Color = st_cmap(stc);

         he(stc,dc).CapSize = 6;
         
        he(stc,dc).LineWidth = lw;
        he(stc,dc).LineStyle = 'none';
        he(stc,dc).DisplayName = MergedData(dc,stc,1).donors_str{:};
        
        
        matforpvalues(dc,stc) = delta2_mu;
        errmatforpvalues(dc,stc) = 2*err_delta2_mu; %twice because err_delta_mu was halved for plotting
        
    end %for dc
end %for

save('lambda2_ratio(48h)_CtrlNormalised.mat',...
    'delta2_lambda2_48h_mat',...
    'le_delta2_lambda2_48h_mat',...
    'ue_delta2_lambda2_48h_mat',...
    'sampletypes');

setsemilogy




% axes properties
%ha.XTick = 1:numel(sampletypes)-1;
ha.XTick = 2:numel(sampletypes);

%ha.XLim = minmax(ha.XTick) + 0.5.*[-1 1];
ha.XLim = [1.2,3.2];%minmax(ha.XTick) + 0.5.*[-1 1];
ha.YLim = [1e-2,1e2];%minmax(ha.XTick) + 0.5.*[-1 1];

ha.XTickLabel = {'VX809-VX770','RPL554'};
ha.XTickLabelRotation = 15;

plot(xlim,ones(2,1),'--','color',[0.5 0.5 0.5],'linewidth',0.4)

% ha.XLabel.String = ' ';
% ha.XLabel.FontSize = lfs;
ha.YLabel.String = '(\lambda^2_{48h} / \lambda^2_{0h}) / (\lambda^2_{48h} / \lambda^2_{0h})_{ctrl}';
% ht = fix_xticklabels(ha);
% set(ht,'FontSize',lfs,...
%     'FontName','Arial',...
%     'Color',darkGrey);

hl = legend(he(end,:));
hl.String = replace(hl.String, {'_L';'_R'},{'Lumacaftor + Ivacaftor';'RPL554'});
hl.Box = 'off';
hl.FontSize = tfs*2;
hl.Units = 'normalized';
hl.String = regexprep(hl.String,'d(?=\d)','s');
% positions
%ha.Position = hapos + [0.03 0 0 0];
%hl.Position(1) = ha.Position(1) - 0.05;
%hl.Position(2) = ha.Position(2) + ha.Position(4) - hl.Position(4) - 0.01;
hl.Location = 'northwest'
ha.YLabel.Units = 'normalized';
ha.YLabel.Position(1) = -0.13;
ha.YLabel.VerticalAlignment = 'baseline';


savename = 'lambda2_ratio_48h_normalised_by_control.pdf';
%print2svg(hf, savename, 1);
%print2pdf(hf, savename, 0);
saveas(gcf,savename)
% 
% disp('using paired t-test')
% for stc = 2:numel(sampletypes)
%     disp(sampletypes{stc})
%     disp('P value \lambda2(48h) / \lambda2(0h) vs control');
%     [~,p] = ttest(matforpvalues(:,stc), 0,'tail','left');
%     disp(p)
% end



matforpvalues
errmatforpvalues


%% %%% get values of CBF(48h) - CBF(00h) - (same on DMSO) and plot
% it's the same things we did on lambda2 to confront CFCF with CFnull, but
% on CBF

bsi = find(boxsizes_vector == 64);
delta2_CBF_48h_mat = [];
e_delta2_CBF_48h_mat = [];
for stc = 2:numel(sampletypes)
    for dc = 1:numel(donors)
        
        % no need to apply the left shift as we are then subtracting one
        % from the other
        CBF_48h = nanmean(MergedData(dc,stc,2).Frequency_Hz{bsi});
        CBF_48h_std = nanstd(MergedData(dc,stc,2).Frequency_Hz{bsi});
        CBF_00h = nanmean(MergedData(dc,stc,1).Frequency_Hz{bsi});
        CBF_00h_std = nanstd(MergedData(dc,stc,1).Frequency_Hz{bsi});
        
        % normalise by 00h and propagate error
        delta_CBF = CBF_48h - CBF_00h;
        err_delta_CBF = sqrt(CBF_48h_std.^2 + CBF_00h_std.^2);
        
        % now same for control (yeah I shouldn;t recalculate it every time
        % but whatever)
        CBF_48h_ctrl = nanmean(MergedData(dc,1,2).Frequency_Hz{bsi});
        CBF_48h_std_ctrl = nanstd(MergedData(dc,1,2).Frequency_Hz{bsi});
        CBF_00h_ctrl = nanmean(MergedData(dc,1,1).Frequency_Hz{bsi});
        CBF_00h_std_ctrl = nanstd(MergedData(dc,1,1).Frequency_Hz{bsi});
        
        % normalise ctrl by 00h and propagate error
        delta_CBF_ctrl = CBF_48h_ctrl - CBF_00h_ctrl;
        err_delta_CBF_ctrl = sqrt(CBF_48h_std_ctrl.^2 + CBF_00h_std_ctrl.^2);
        
        % normalise by control
        delta2_CBF = delta_CBF - delta_CBF_ctrl;
        err_delta2_CBF = sqrt(err_delta_CBF.^2 + err_delta_CBF_ctrl.^2);
        
        % things are already in the right units
                 
        
        delta2_CBF_48h_mat(dc,stc) = delta2_CBF; %y to plot
        e_delta2_CBF_48h_mat(dc,stc) = err_delta2_CBF; %error        
        
        
    end %for dc
end %for


save('normDeltaCBF(48h).mat',...
    'delta2_CBF_48h_mat',...
    'e_delta2_CBF_48h_mat',...
    'sampletypes');

% now plot

markarr = '^odsv<>';
linest = {'-','-','--','--','-.'};
darkGrey = [0.15 0.15 0.15];
st_cmap = [27,158,119
    217,95,2
    152,67,1
    117,112,179
    231,41,138
    102,166,30 ] ./ 255;
st_cmap = vertcat(st_cmap, darkGrey);
x_stagger = [-0.09 0 0.09];
lfs = 10;   % labels fontsize
tfs = 8;    % ticks fontsize
lw = 1.2;

%[hf, ha] = figure_template;

st_cmap=[1 0 0
    1 0 0
    0 0 1];

[hf, ha] = figure_template;
%hf = figure;
%ha = axes;



clear he;

% sampletype on x axes, 
%sampletype_order = [5 4 6 1 2 3];
sampletype_order = [1 2 3];
for stc = 2:numel(sampletypes)
    for dc = 1:numel(donors)
        
        sto = sampletype_order(stc)
        he(stc,dc) = errorbar(stc + x_stagger(dc),...
            delta2_CBF_48h_mat(dc,sto),...
            e_delta2_CBF_48h_mat(dc,sto));

        he(stc,dc).Marker = markarr(dc);
        he(stc,dc).MarkerFaceColor = 'w';
        he(stc,dc).Color = st_cmap(sto,:);
        he(stc,dc).CapSize = 6;
        he(stc,dc).LineWidth = lw;
        he(stc,dc).LineStyle = 'none';
        he(stc,dc).DisplayName = MergedData(dc,sto,1).donors_str{:};
        
    end %for dc
end %for


% axes properties
ha.XTick = 2:numel(sampletypes);
%ha.XLim = minmax(ha.XTick) + 0.5.*[-1 1];
ha.XLim = [1.5 3.5];

%ha.XTickLabel = {MergedData(1,sampletype_order,1).display_name};
ha.XTickLabel = {'VX809-VX770','RPL554'};

ha.XTickLabelRotation = 15;
% ha.YLim = [-10 6];
ha.YLim = [-30 30];

plot(xlim,zeros(2,1),'--','color',[0.5 0.5 0.5],'linewidth',0.4)

% ha.XLabel.String = ' ';
% ha.XLabel.FontSize = lfs;
ha.YLabel.String = '\Delta_{CBF}(48h) - \Delta_{CBF}^{ctrl}(0h)';
% ht = fix_xticklabels(ha);
% set(ht,'FontSize',lfs,...
%     'FontName','Arial',...
%     'Color',darkGrey);

hl = legend(he(end,:));
hl.Box = 'off';
hl.FontSize = tfs;
hl.Units = 'normalized';
hl.String = regexprep(hl.String,'d(?=\d)','s');
% positions
%ha.Position = hapos + [0.03 0 0 0];
hl.Position(1) = ha.Position(1) - 0.05;
hl.Position(2) = ha.Position(2) + 0.01;
hl.Location='northeast'
ha.YLabel.Units = 'normalized';
ha.YLabel.Position(1) = -0.13;
ha.YLabel.VerticalAlignment = 'baseline';


savename = 'normDeltaCBF(48h)_errorbar';
%print2svg(hf, savename, 1);
%print2pdf(hf, savename, 0);

% disp('using paired t-test')
% for stc = 1:numel(sampletypes)-1
%     disp(sampletypes{stc})
%     disp('P value CBF(48h) - CBF(0h) vs control');
%     [~,p] = ttest(delta2_CBF_48h_mat(:,stc), 0,'tail','right');
%     disp(p)
% end
% 
% 

%%

%% %%% get values of CBF(48h) - CBF(00h) - (same on DMSO) and plot
% same as before but with plot on donors

bsi = find(boxsizes_vector == 64);
delta2_CBF_48h_mat = [];
e_delta2_CBF_48h_mat = [];
for stc = 2:numel(sampletypes)
    for dc = 1:numel(donors)
        
        % no need to apply the left shift as we are then subtracting one
        % from the other
        CBF_48h = nanmedian(MergedData(dc,stc,2).Frequency_Hz{bsi});
        CBF_48h_std = nanstd(MergedData(dc,stc,2).Frequency_Hz{bsi});
        CBF_00h = nanmedian(MergedData(dc,stc,1).Frequency_Hz{bsi});
        CBF_00h_std = nanstd(MergedData(dc,stc,1).Frequency_Hz{bsi});
        
        % normalise by 00h and propagate error
        delta_CBF = CBF_48h - CBF_00h;
        err_delta_CBF = sqrt(CBF_48h_std.^2 + CBF_00h_std.^2);
        
        % now same for control (yeah I shouldn;t recalculate it every time
        % but whatever)
        CBF_48h_ctrl = nanmedian(MergedData(dc,1,2).Frequency_Hz{bsi});
        CBF_48h_std_ctrl = nanstd(MergedData(dc,1,2).Frequency_Hz{bsi});
        CBF_00h_ctrl = nanmedian(MergedData(dc,1,1).Frequency_Hz{bsi});
        CBF_00h_std_ctrl = nanstd(MergedData(dc,1,1).Frequency_Hz{bsi});
        
        % normalise ctrl by 00h and propagate error
        delta_CBF_ctrl = CBF_48h_ctrl - CBF_00h_ctrl;
        err_delta_CBF_ctrl = sqrt(CBF_48h_std_ctrl.^2 + CBF_00h_std_ctrl.^2);
        
        % normalise by control
        delta2_CBF = delta_CBF - delta_CBF_ctrl;
        err_delta2_CBF = sqrt(err_delta_CBF.^2 + err_delta_CBF_ctrl.^2);
        
        % things are already in the right units
                 
        
        delta2_CBF_48h_mat(dc,stc) = delta2_CBF; %y to plot
        e_delta2_CBF_48h_mat(dc,stc) = err_delta2_CBF; %error        
        
        
    end %for dc
end %for


save('normDeltaCBF(48h).mat',...
    'delta2_CBF_48h_mat',...
    'e_delta2_CBF_48h_mat',...
    'sampletypes');

% now plot

markarr = '^odsv<>';
linest = {'-','-','--','--','-.'};
darkGrey = [0.15 0.15 0.15];
st_cmap = [27,158,119
    217,95,2
    152,67,1
    117,112,179
    231,41,138
    102,166,30 ] ./ 255;
st_cmap = vertcat(st_cmap, darkGrey);
x_stagger = [-0.09 0 0.09];
lfs = 10;   % labels fontsize
tfs = 8;    % ticks fontsize
lw = 1.2;

%[hf, ha] = figure_template;

st_cmap=[1 0 0
    1 0 0
    0 0 1];

%[hf, ha] = figure_template;
hf = figure;
ha = axes;



clear he;

% sampletype on x axes, 
%sampletype_order = [5 4 6 1 2 3];
sampletype_order = [1 2 3];
for stc = 2:numel(sampletypes)
    for dc = 1:numel(donors)
        
        sto = sampletype_order(stc)
        he(stc,dc) = errorbar(dc + x_stagger(stc),...
            delta2_CBF_48h_mat(dc,sto),...
            e_delta2_CBF_48h_mat(dc,sto));
        hold on;

        he(stc,dc).Marker = markarr(dc);
        he(stc,dc).MarkerFaceColor = 'w';
        he(stc,dc).Color = st_cmap(sto,:);
        he(stc,dc).CapSize = 6;
        he(stc,dc).LineWidth = lw;
        he(stc,dc).LineStyle = 'none';
        he(stc,dc).DisplayName = MergedData(dc,sto,1).sampletype_str;
        
    end %for dc
end %for


% axes properties
ha.XTick = 1:numel(sampletypes);
%ha.XLim = minmax(ha.XTick) + 0.5.*[-1 1];
ha.XLim = [0.5 3.5];

%ha.XTickLabel = {MergedData(1,sampletype_order,1).display_name};
ha.XTickLabel = {'CF1','CF2','CF3'};

ha.XTickLabelRotation = 15;
% ha.YLim = [-10 6];
ha.YLim = [-15 15];

plot(xlim,zeros(2,1),'--','color',[0.5 0.5 0.5],'linewidth',0.4)

% ha.XLabel.String = ' ';
% ha.XLabel.FontSize = lfs;
ha.YLabel.String = '\Delta_{CBF}(48h) - \Delta_{CBF}^{ctrl}(0h)';
% ht = fix_xticklabels(ha);
% set(ht,'FontSize',lfs,...
%     'FontName','Arial',...
%     'Color',darkGrey);

hl = legend(he(2:end,end));
hl.String = replace(hl.String, {'_L';'_R'},{'Lumacaftor + Ivacaftor';'RPL554'});
%hl = legend('VX809-VX770','RPL554');
hl.Box = 'on';
hl.FontSize = tfs*1.7;
hl.Units = 'normalized';
%hl.String = regexprep(hl.String,'d(?=\d)','s');
% positions
%ha.Position = hapos + [0.03 0 0 0];
hl.Position(1) = ha.Position(1) - 0.05;
hl.Position(2) = ha.Position(2) + 0.01;
hl.Location='northeast'
ha.YLabel.Units = 'normalized';
ha.YLabel.Position(1) = -0.13;
ha.YLabel.VerticalAlignment = 'baseline';
ha.XLabel.String = 'Donor';

savename = 'normDeltaCBF(48h)_errorbar_donors.pdf';
saveas(gcf,savename)
%print2svg(hf, savename, 1);
%print2pdf(hf, savename, 0);

% disp('using paired t-test')
% for stc = 1:numel(sampletypes)-1
%     disp(sampletypes{stc})
%     disp('P value CBF(48h) - CBF(0h) vs control');
%     [~,p] = ttest(delta2_CBF_48h_mat(:,stc), 0,'tail','right');
%     disp(p)
% end
% 
% 


%% CBF boxplots SEPARATED DONORS vs control

darkB = [0 0.4470 0.7410];
darkO = [0.8500 0.3250 0.0980];
darkG = [0.13 0.53 0.13];

lightB = darkB + 0.25;
lightO = darkO + [0.15 0.23 0.2];
darkGrey = [0.15 0.15 0.15];
lightGrey = [0.65 0.65 0.65];
cmap = [27,158,119 ;
    217,95,2 ;
    117,112,179 ] ./ 255;

st_cmap = vertcat(darkB, darkG, darkO, darkGrey) ;

st_cmap =[ 0,0,0;
    1,0,0;
    0,0,1;]
%%close all

lfs = 10;   % labels fontsize
tfs = 8;    % ticks fontsize
lw = 1.2;

bsi = find(boxsizes_vector == 64);

% main for loop on donors

for dc = 1:numel(donors)
    
    
    % prepare data for boxplot
    
    CBF = [];   % CBF
    ST = [];    % 1st grouping variable
    TP = [];    % 2nd grouping variable
    
    for stc = numel(sampletypes):-1:1
        
        for tpc = 1:numel(timepoints)
            
            dummy = MergedData(dc, stc,tpc).Frequency_Hz{bsi};
            
            ST = vertcat(ST, repmat(sampletypes(stc), size(dummy)));
            TP = vertcat(TP, repmat(timepoints(tpc),size(dummy)));
            CBF = vertcat(CBF, dummy);
            
        end %tpc
    end %for stc
    
    
    time = 1:2;%24 .* (0:3);
    % now plot
 
        
        idx_ctrl = contains(ST, '_D' ); %select ctrl and the sampletype we want to plot
         
        hf = figure;
        hf.Units = 'centimeters';
        hf.Position([3 4]) = [8 6];
        
        ha = axes;
        ha.Box = 'on';
        ha.NextPlot = 'add';
        ha.FontSize = tfs;
        ha.FontName = 'Arial';
        
        pos_vec = reshape(repmat(time, 2, 1) + [-0.1 0.1]',[],1);
        pos_vec= [0.8,1,1.2,1.8,2,2.2];
        % ctrl first
 %       hb_ctrl = boxplot(CBF(idx_ctrl),{TP(idx_ctrl), ST(idx_ctrl)},...
       hb_ctrl = boxplot(CBF(idx_ctrl),strcat(TP(idx_ctrl), ST(idx_ctrl)),...
        'ColorGroup',ST(idx_ctrl),'Positions',pos_vec([1,4]),...
            'Symbol','+',...
            'OutlierSize',3,...
            'Colors',darkGrey,...
            'Widths',0.15);
        ha.XTickLabel = {''};
        
     for stc = 2:numel(sampletypes)      
       
         idx_drug = contains(ST, sampletypes{stc} ); %select ctrl and the sampletype we want to plot

        
        % then drug
        hb_drug = boxplot(CBF(idx_drug),{TP(idx_drug), ST(idx_drug)},...
            'ColorGroup',ST(idx_drug),'Positions',pos_vec([stc,stc+3]),...
            'Symbol','+',...
            'OutlierSize',3,...
            'Colors',st_cmap(stc,:),...
            'Widths',0.15);
        ha.XTickLabel = {''};
        
        
       end %for stc
        
        drawnow
        
        ha.XTickLabel = {''};
        ha.XTick = time;
        ha.XTickLabel = timepoints;
        ha.FontSize = tfs;
        
        ha.XLim = [0.5 numel(timepoints) + 0.5];
        ha.YLim = [0 15];
        
        %     ha.Position = [0.13 0.18 0.85 0.79];
        
        ha.XLabel.String = 'Time';
        ha.XLabel.FontSize = lfs;
        
        ha.YLabel.String = 'CBF, [Hz]';
        ha.YLabel.FontSize = lfs;
        
        % fix title
 %       ha.Title.String = strjoin({donors{dc},sampletypes{[1 stc]} },' ');
        ha.Title.String = {donors{dc}};

 %       ha.Title.String = replace(ha.Title.String, {'LumacaFTOR';'Orkambi';'T1a'},{'Lumacaftor';'Lumacaftor + Ivacaftor';'T\alpha1'});
        ha.Title.FontWeight = 'normal';
        ha.Title.FontSize = lfs;
        
        box_handles = findall(ha,'Tag','Box');
        
%         set(findall(ha,'tag','Box'),'LineWidth',lw);
%         set(findall(ha,'tag','Median'),'LineWidth',lw);
%         set(findall(ha,'tag','Lower Whisker'),'LineWidth',lw);
%         set(findall(ha,'tag','Upper Whisker'),'LineWidth',lw);
%         set(findall(ha,'tag','Lower Adjacent Value'),'LineWidth',lw);
%         set(findall(ha,'tag','Upper Adjacent Value'),'LineWidth',lw);
          alllines=findall(gca,'Tag','Box');
    hLegend = legend(alllines(1:2:end), {'RPL554','Lumacaftor+Ivacaftor','DMSO'});
          %             hleg = legend(box_handles(1),sampletypes(1)' );
        %     hleg.Box = 'off';
        %     hleg.FontSize = lfs;
        %     hleg.Interpreter = 'none';
        %
        %     hleg.Position([1 2]) = [ha.Position(1) sum(ha.Position([2 4])) - 1.15*hleg.Position(4)];
        
        savename = strjoin({'CBF_CtrlBoxplot',...
            [donors{dc}],'.pdf'});
        disp(savename)
 %       savename = fullfile(figures_folder,savename);
        saveas(gcf,savename);
%        print2svg(hf, savename, 1, 1);
        
    
    
end % for dc


