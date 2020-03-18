% D:\Dropbox\Universita\PhD\Papers\DDM_and_profiles_paper\Figures\fig4_CF_drugs_newboxes


cch

cd('D:\Dropbox\Universita\PhD\Papers\DDM_and_profiles_paper\Figures\fig4_CF_drugs_newboxes');

%% plot frequencies over time, one per donor

cch
clearvars -except dc
% load('D:\Data\Cilia\Analysis\2017_08_04--07\AccumData_from_func_usedfordraft.mat',...
load('D:\Data\Cilia\Analysis\2017_08_04--07\AccumData_from_func_newboxes.mat',...
    'SampleType',...
    'figures_folder',..._nosparse
    'sampletypes',...
    'timepoints',...
    'donors',...
    'inserts',...
    'positions',...
    'boxsizes_vector',...
    'q_limits_1oum');



for stc = 4:-1:1
    for tpc = 4:-1:1
        
        MergedData_alldonors(stc, tpc) = merge_SampleType_data(SampleType,...
                stc, tpc, donors, '', '');
        
        % separate donors: we have 3 donors, 4 drugs, 4 timepoints
        for dc = 3:-1:1
            
            MergedData(dc,stc,tpc) = merge_SampleType_data(SampleType,...
                stc, tpc, donors{dc}, '', '');
            
        end %for dc
        
    end %for tpc
end %for stc

figures_folder = 'D:\Dropbox\Universita\PhD\Papers\DDM_and_profiles_paper\Figures\fig4_CF_drugs_newboxes';


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
    ha.Title.String = regexprep(ha.Title.String,'d(?=\d)','s');
    
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
    %print2svg(hf, savename, 1, 1)
    
end

%% save CBF values in table

% 4 tables horizonatally (one per treatment)
% 3 rows (one per timepoint)
alphabet = 'A':'Z';
alphabet2 = [' ',alphabet];


clear xlscolumnsname
for i = 1:600
    fci = ceil(i/26);
    sci = mod(i-1,26)+1;
    xlscolumnsname{i,1} = [alphabet2(fci),alphabet(sci)];
end

%%
% one table per timepoint
for tpc = 1:numel(timepoints)
    for stc = 1:numel(sampletypes)
        
        T = table(...
            regexprep(donors,'d(?=\d)','s'),...
            arrayfun(@(i) mean(MergedData(i,stc,tpc).Frequency_Hz{bsi}), 1:3)', ...
            arrayfun(@(i) std(MergedData(i,stc,tpc).Frequency_Hz{bsi}), 1:3)', ...
            arrayfun(@(i) prctile(MergedData(i,stc,tpc).Frequency_Hz{bsi},5), 1:3)', ...
            arrayfun(@(i) prctile(MergedData(i,stc,tpc).Frequency_Hz{bsi},95), 1:3)', ...
            repmat({'Hz'},3,1),...
            {MergedData(:,stc,tpc).sampletype_str}',...
            {MergedData(:,stc,tpc).timepoint_str}',...
            'VariableNames',{'subject','average','stdev','prctile_5th','prctile_95th','units','sampletype','timepoint'}...
            );
        
        disp(T)
        nletters = ceil(((stc-1) * (size(T,2)+2) + 1) / 26);
        range_letter = xlscolumnsname{(stc-1) * (size(T,2)+2) + 1,1};
        range_number = num2str((tpc-1) * (numel(donors)+3) +1 );

%         writetable(T,'fig4_CBF_table.xlsx', 'Sheet', 1, 'Range',[range_letter,range_number],...
%                 'WriteVariableNames',true)
    end
end

% winopen('fig4_CBF_table.xlsx')



%% now figure of histogram of control + one drug at a time

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
binedges = 0:0.5:30;
%     binedges = 0:0.4:15;

bsi = find(boxsizes_vector == 64);

%%close all


for dc = 1:numel(donors)
    for stc = 1:numel(sampletypes)-1
        for tpc = 1:numel(timepoints)
            
            i = sub2ind(size(MergedData), dc,stc,tpc);
            
            hf = figure;
            hf.Units = 'centimeters';
            hf.Position([3 4]) = [8 6];
            
            ha = axes;
            ha.Box = 'on';
            ha.NextPlot = 'add';
            ha.FontSize = tfs;
            ha.FontName = 'Arial';
            
            
            % drug
            histogram(MergedData(i).Frequency_Hz{bsi}, binedges,...
                'Normalization','Probability',...
                'FaceColor', st_cmap(stc,:),...
                'FaceAlpha', 1,...
                'EdgeColor', 'none');
            
            % control
            histogram(MergedData(dc,4,tpc).Frequency_Hz{bsi}, binedges,...
                'Normalization','Probability',...
                'DisplayStyle','Stairs',...
                'EdgeColor', darkGrey,...
                'LineStyle', '-',...
                'LineWidth',1.5);
            
            ha.XLabel.String = 'CBF, [Hz]';
            ha.XLabel.FontSize = lfs;
            
            ha.YLabel.String = 'Counts, normalised';
            ha.YLabel.FontSize = lfs;
            
            % fix title
            ha.Title.String = strjoin({MergedData(i).donors_str{:},MergedData(i).timepoint_str},' ');
            ha.Title.FontWeight = 'normal';
            ha.Title.FontSize = lfs;
            ha.Title.String = regexprep(ha.Title.String,'d(?=\d)','s');

            ha.YLim = [0 0.45];
            ha.XLim = [0 20];
            
            hleg = legend({MergedData(i).sampletype_str, MergedData(dc,4,tpc).sampletype_str});
            hleg.String = replace(hleg.String, {'LumacaFTOR';'Orkambi';'T1a'},{'Lumacaftor';'Lumacaftor + Ivacaftor';'T\alpha1'});
            hleg.Box = 'off';
            hleg.FontSize = lfs;
%             hleg.Interpreter = 'none';
            
            hleg.Position([1 2]) = [sum(ha.Position([1 3])) - hleg.Position(3) sum(ha.Position([2 4])) - 1.15*hleg.Position(4)];
            
            
            savename = strjoin({'CBF_CtrlHistogram',...
                MergedData(i).donors_str{:},...
                MergedData(i).sampletype_str,...
                MergedData(i).timepoint_str},'_');
            savename = fullfile(figures_folder,savename);
            %print2svg(hf, savename, 1, 1);
            
            
        end %for tpc
    end %for stc
end %for dc


%% now figure of histograms ALL DONORS 

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


for i = 1:numel(MergedData_alldonors) 
    
    [stc,~] = ind2sub(  size(MergedData_alldonors),i);

    
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
    histogram(MergedData_alldonors(i).Frequency_Hz{bsi}, binedges,...
        'Normalization','Probability',...
        'FaceColor', st_cmap(stc,:),...
        'FaceAlpha', 1,...
        'EdgeColor', 'none');
    
    ha.XLabel.String = 'CBF, [Hz]';
    ha.XLabel.FontSize = lfs;
    
    ha.YLabel.String = 'Counts, normalised';
    ha.YLabel.FontSize = lfs;
    
    % fix title
    ha.Title.String = strjoin({MergedData_alldonors(i).donors_str{:},...
        MergedData_alldonors(i).timepoint_str},' ');
    ha.Title.FontWeight = 'normal';
    ha.Title.FontSize = lfs;
    ha.Title.String = regexprep(ha.Title.String,'d(?=\d)','s');

    ha.YLim = [0 0.3];
    ha.XLim = [0 20];

    hleg = legend({MergedData_alldonors(i).sampletype_str});
    hleg.String = replace(hleg.String,  {'LumacaFTOR';'Orkambi';'T1a'},{'Lumacaftor';'Lumacaftor + Ivacaftor';'T\alpha1'});
    hleg.Box = 'off';
    hleg.FontSize = lfs;
%     hleg.Interpreter = 'none';
    
    hleg.Position([1 2]) = [sum(ha.Position([1 3])) - hleg.Position(3) sum(ha.Position([2 4])) - 1.3*hleg.Position(4)];
    
    
    savename = strjoin({'CBF_histogram',...
        [MergedData_alldonors(i).donors_str{:}],...
        MergedData_alldonors(i).sampletype_str,...
        MergedData_alldonors(i).timepoint_str},'_');
    savename = fullfile(figures_folder,savename);
%print2svg(hf, savename, 1, 1)
    
end


%% now figure of histogram  ALL DONORS of control + one drug at a time

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
binedges = 0:0.5:30;
%     binedges = 0:0.4:15;

bsi = find(boxsizes_vector == 64);

%%close all


for stc = 1:numel(sampletypes)-1
    for tpc = 1:numel(timepoints)
        
        i = sub2ind(size(MergedData_alldonors), stc,tpc);
        
        hf = figure;
        hf.Units = 'centimeters';
        hf.Position([3 4]) = [8 6];
        
        ha = axes;
        ha.Box = 'on';
        ha.NextPlot = 'add';
        ha.FontSize = tfs;
        ha.FontName = 'Arial';
        
        
        % drug
        histogram(MergedData_alldonors(i).Frequency_Hz{bsi}, binedges,...
            'Normalization','Probability',...
            'FaceColor', st_cmap(stc,:),...
            'FaceAlpha', 1,...
            'EdgeColor', 'none');
        
        % control
        histogram(MergedData_alldonors(4,tpc).Frequency_Hz{bsi}, binedges,...
            'Normalization','Probability',...
            'DisplayStyle','Stairs',...
            'EdgeColor', darkGrey,...
            'LineStyle', '-',...
            'LineWidth',1.5);
        
        ha.XLabel.String = 'CBF, [Hz]';
        ha.XLabel.FontSize = lfs;
        
        ha.YLabel.String = 'Counts, normalised';
        ha.YLabel.FontSize = lfs;
        
        % fix title
        ha.Title.String = strjoin({MergedData_alldonors(i).donors_str{:},MergedData_alldonors(i).timepoint_str},' ');
        ha.Title.FontWeight = 'normal';
        ha.Title.FontSize = lfs;
        ha.Title.String = regexprep(ha.Title.String,'d(?=\d)','s');

        ha.YLim = [0 0.33];
        ha.XLim = [0 20];
        
        hleg = legend({MergedData_alldonors(i).sampletype_str, MergedData_alldonors(4,tpc).sampletype_str});
        hleg.String = replace(hleg.String,  {'LumacaFTOR';'Orkambi';'T1a'},{'Lumacaftor';'Lumacaftor + Ivacaftor';'T\alpha1'});
        hleg.Box = 'off';
        hleg.FontSize = lfs;
%         hleg.Interpreter = 'none';
        
        hleg.Position([1 2]) = [sum(ha.Position([1 3])) - hleg.Position(3) sum(ha.Position([2 4])) - 1.15*hleg.Position(4)];
        
        
        savename = strjoin({'CBF_CtrlHistogram',...
            [MergedData_alldonors(i).donors_str{:}],...
            MergedData_alldonors(i).sampletype_str,...
            MergedData_alldonors(i).timepoint_str},'_');
        savename = fullfile(figures_folder,savename);
        %print2svg(hf, savename, 1, 1);
        
        
    end %for tpc
end %for stc

%% CBF boxplots ALL DONORS vs control


close all

flag_no72h = 0;

lfs = 10;   % labels fontsize
tfs = 8;    % ticks fontsize
lw = 1.2;

bsi = find(boxsizes_vector == 64);

% prepare data for boxplot

CBF = [];   % CBF
ST = [];    % 1st grouping variable
TP = [];    % 2nd grouping variable

for stc = numel(sampletypes):-1:1
    
    for tpc = 1:numel(timepoints)
        
        dummy = MergedData_alldonors(stc,tpc).Frequency_Hz{bsi};
        
        ST = vertcat(ST, repmat(sampletypes(stc), size(dummy)));
        TP = vertcat(TP, repmat(timepoints(tpc),size(dummy)));
        CBF = vertcat(CBF, dummy);
        
    end %tpc
end %for stc
    
if flag_no72h
    time = 1:3;
else
    time = 1:4;%24 .* (0:3);
end
% now plot
for stc = 1:numel(sampletypes) - 1
    
    idx_ctrl = contains(ST, 'Control' ); %select ctrl and the sampletype we want to plot
    idx_drug = contains(ST, sampletypes{stc} ); %select ctrl and the sampletype we want to plot
    
    % take care of timing
    idx_ctrl = idx_ctrl & contains(TP,timepoints(time));
    idx_drug = idx_drug & contains(TP,timepoints(time));
    
    hf = figure;
    hf.Units = 'centimeters';
    hf.Position([3 4]) = [8 6];
    
    ha = axes;
    ha.Box = 'on';
    ha.NextPlot = 'add';
    ha.FontSize = tfs;
    ha.FontName = 'Arial';
    
    pos_vec = reshape(repmat(time, 2, 1) + [-0.1 0.1]',[],1);

    % ctrl first
    hb_ctrl = boxplot(CBF(idx_ctrl),{TP(idx_ctrl), ST(idx_ctrl)},...
        'ColorGroup',ST(idx_ctrl),'Positions',pos_vec(1:2:end),...
        'Symbol','+',...
        'OutlierSize',3,...
        'Colors',darkGrey,...
        'Widths',0.15);
    ha.XTickLabel = {''};

    % then drug
    hb_drug = boxplot(CBF(idx_drug),{TP(idx_drug), ST(idx_drug)},...
        'ColorGroup',ST(idx_drug),'Positions',pos_vec(2:2:end),...
        'Symbol','+',...
        'OutlierSize',3,...
        'Colors',st_cmap(stc,:),...
        'Widths',0.15);
        ha.XTickLabel = {''};

    drawnow
    
    ha.XTickLabel = {''};
    ha.XTick = time;
    ha.XTickLabel = timepoints(time);
    ha.FontSize = tfs;
    
    ha.XLim = [0.5 numel(timepoints(time)) + 0.5];
    ha.YLim = [0 30];
    
    %     ha.Position = [0.13 0.18 0.85 0.79];
    
    ha.XLabel.String = 'Time';
    ha.XLabel.FontSize = lfs;
    
    ha.YLabel.String = 'CBF, [Hz]';
    ha.YLabel.FontSize = lfs;
    
    % fix title
    ha.Title.String = strjoin({'d1','d2','d3',sampletypes{[4 stc]} },' ');
    ha.Title.String = replace(ha.Title.String, {'LumacaFTOR';'Orkambi';'T1a'},{'Lumacaftor';'Lumacaftor + Ivacaftor';'T\alpha1'});
    ha.Title.FontWeight = 'normal';
    ha.Title.FontSize = lfs;
    ha.Title.String = regexprep(ha.Title.String,'d(?=\d)','s');

    box_handles = findall(ha,'Tag','Box');
    
%     set(findall(ha,'tag','Box'),'LineWidth',lw);
%     set(findall(ha,'tag','Median'),'LineWidth',lw);
%     set(findall(ha,'tag','Lower Whisker'),'LineWidth',lw);
%     set(findall(ha,'tag','Upper Whisker'),'LineWidth',lw);
%     set(findall(ha,'tag','Lower Adjacent Value'),'LineWidth',lw);
%     set(findall(ha,'tag','Upper Adjacent Value'),'LineWidth',lw);

    %     hleg = legend(box_handles([7 3]),sampletypes([4 stc])' );
    %     hleg.Box = 'off';
    %     hleg.FontSize = lfs;
    %     hleg.Interpreter = 'none';
    %
    %     hleg.Position([1 2]) = [ha.Position(1) sum(ha.Position([2 4])) - 1.15*hleg.Position(4)];
    
    savename = strjoin({'CBF_CtrlBoxplot',...
        [donors{:}],...
        sampletypes{stc}},'_');
    if flag_no72h
        savename = [savename,'_no72h'];
    end
    savename = fullfile(figures_folder,savename);
    %print2svg(hf, savename, 1, 1);    
end %for stc

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

%%close all

flag_no72h = 1;


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
    
    
    if flag_no72h
        time = 1:3;
    else
        time = 1:4;%24 .* (0:3);
    end
    % now plot
    for stc = 1:numel(sampletypes) - 1
        
        idx_ctrl = contains(ST, 'Control' ); %select ctrl and the sampletype we want to plot
        idx_drug = contains(ST, sampletypes{stc} ); %select ctrl and the sampletype we want to plot
        
        % take care of timing
        idx_ctrl = idx_ctrl & contains(TP,timepoints(time));
        idx_drug = idx_drug & contains(TP,timepoints(time));
        
        hf = figure;
        hf.Units = 'centimeters';
        hf.Position([3 4]) = [8 6];
        
        ha = axes;
        ha.Box = 'on';
        ha.NextPlot = 'add';
        ha.FontSize = tfs;
        ha.FontName = 'Arial';
        
        pos_vec = reshape(repmat(time, 2, 1) + [-0.1 0.1]',[],1);
        
        % ctrl first
        hb_ctrl = boxplot(CBF(idx_ctrl),{TP(idx_ctrl), ST(idx_ctrl)},...
            'ColorGroup',ST(idx_ctrl),'Positions',pos_vec(1:2:end),...
            'Symbol','+',...
            'OutlierSize',3,...
            'Colors',darkGrey,...
            'Widths',0.15);
        ha.XTickLabel = {''};
        
        % then drug
        hb_drug = boxplot(CBF(idx_drug),{TP(idx_drug), ST(idx_drug)},...
            'ColorGroup',ST(idx_drug),'Positions',pos_vec(2:2:end),...
            'Symbol','+',...
            'OutlierSize',3,...
            'Colors',st_cmap(stc,:),...
            'Widths',0.15);
        ha.XTickLabel = {''};
        
        drawnow
        
        ha.XTickLabel = {''};
        ha.XTick = time;
        ha.XTickLabel = timepoints(time);
        ha.FontSize = tfs;
        
        ha.XLim = [0.5 numel(timepoints(time)) + 0.5];
        ha.YLim = [0 30];
        
        %     ha.Position = [0.13 0.18 0.85 0.79];
        
        ha.XLabel.String = 'Time';
        ha.XLabel.FontSize = lfs;
        
        ha.YLabel.String = 'CBF, [Hz]';
        ha.YLabel.FontSize = lfs;
        
        % fix title
        ha.Title.String = strjoin({donors{dc},sampletypes{[4 stc]} },' ');
        ha.Title.String = replace(ha.Title.String, {'LumacaFTOR';'Orkambi';'T1a'},{'Lumacaftor';'Lumacaftor + Ivacaftor';'T\alpha1'});
        ha.Title.FontWeight = 'normal';
        ha.Title.FontSize = lfs;
        ha.Title.String = regexprep(ha.Title.String,'d(?=\d)','s');
        
        box_handles = findall(ha,'Tag','Box');
        
%         set(findall(ha,'tag','Box'),'LineWidth',lw);
%         set(findall(ha,'tag','Median'),'LineWidth',lw);
%         set(findall(ha,'tag','Lower Whisker'),'LineWidth',lw);
%         set(findall(ha,'tag','Upper Whisker'),'LineWidth',lw);
%         set(findall(ha,'tag','Lower Adjacent Value'),'LineWidth',lw);
%         set(findall(ha,'tag','Upper Adjacent Value'),'LineWidth',lw);
        
        %     hleg = legend(box_handles([7 3]),sampletypes([4 stc])' );
        %     hleg.Box = 'off';
        %     hleg.FontSize = lfs;
        %     hleg.Interpreter = 'none';
        %
        %     hleg.Position([1 2]) = [ha.Position(1) sum(ha.Position([2 4])) - 1.15*hleg.Position(4)];
        
        savename = strjoin({'CBF_CtrlBoxplot',...
            [donors{dc}],...
            sampletypes{stc}},'_');
        if flag_no72h
            savename = [savename,'_no72h'];
        end
        disp(savename)
        savename = fullfile(figures_folder,savename);
        %print2svg(hf, savename, 1, 1);
        
    end %for stc
    
    
end % for dc

%% CBF boxplots SEPARATED DONORS vs control, all treatments in same axes

close all

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

flag_no72h = 1;


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
    
    for stc = 1:numel(sampletypes)
        
        for tpc = 1:numel(timepoints)
            
            dummy = MergedData(dc, stc,tpc).Frequency_Hz{bsi};
            
            ST = vertcat(ST, repmat(sampletypes(stc), size(dummy)));
            TP = vertcat(TP, repmat(timepoints(tpc),size(dummy)));
            CBF = vertcat(CBF, dummy);
            
        end %tpc
    end %for stc
    
    
    if flag_no72h
        time = 1:3;
    else
        time = 1:4;%24 .* (0:3);
    end
    
    idx_time = contains(TP,timepoints(time));
    
    hf = figure;
    hf.Units = 'centimeters';
    hf.Position([3 4]) = [8 6];
    
    ha = axes;
    ha.Box = 'on';
    ha.NextPlot = 'add';
    ha.FontSize = tfs;
    ha.FontName = 'Arial';
    
    pos_vec = reshape(repmat(time, 4, 1) + linspace(-0.2, 0.2, 4)',[],1);
    
    
    % then drug
    hb_drug = boxplot(CBF(idx_time),{TP(idx_time), ST(idx_time)},...
        'ColorGroup',ST(idx_time),'Positions',pos_vec,...
        'Symbol','+',...
        'OutlierSize',3,...
        'Colors',st_cmap,...
        'Widths',0.08);
    ha.XTickLabel = {''};
    
    drawnow
    
    ha.XTickLabel = {''};
    ha.XTick = time;
    ha.XTickLabel = timepoints(time);
    ha.FontSize = tfs;
    
    ha.XLim = [0.5 numel(timepoints(time)) + 0.5];
    ha.YLim = [0 30];
    
    
    ha.XLabel.String = 'Time';
    ha.XLabel.FontSize = lfs;
    
    ha.YLabel.String = 'CBF, [Hz]';
    ha.YLabel.FontSize = lfs;
    
    % fix title
    ha.Title.String = strjoin({donors{dc}},' ');
    ha.Title.String = replace(ha.Title.String, {'LumacaFTOR';'Orkambi';'T1a'},{'VX809';'VX809 + VX770';'T\alpha1'});
    ha.Title.FontWeight = 'normal';
    ha.Title.FontSize = lfs;
    ha.Title.String = regexprep(ha.Title.String,'d(?=\d)','s');
    
    ha.Position = [0.1591 0.1546 0.82 0.79];

    
    box_handles = findall(ha,'Tag','Box');
    
%             set(findall(ha,'tag','Box'),'LineWidth',lw);
%             set(findall(ha,'tag','Median'),'LineWidth',lw);
%             set(findall(ha,'tag','Lower Whisker'),'LineWidth',lw);
%             set(findall(ha,'tag','Upper Whisker'),'LineWidth',lw);
%             set(findall(ha,'tag','Lower Adjacent Value'),'LineWidth',lw);
%             set(findall(ha,'tag','Upper Adjacent Value'),'LineWidth',lw);
    
    hleg = legend(hb_drug(3,1:4), unique(ST(idx_time),'stable') );
    hleg.Box = 'off';
    hleg.FontSize = lfs;
    hleg.String = replace(hleg.String, {'LumacaFTOR';'Orkambi';'T1a'},{'VX809';'VX809 + VX770';'T\alpha1'});
    hleg.Position([1 2]) = [ha.Position(1) sum(ha.Position([2 4])) - hleg.Position(4)];
    
    savename = strjoin({'CBF_Boxplot',...
        [donors{dc}],...
        'alltreatments'},'_');
    if flag_no72h
        savename = [savename,'_no72h'];
    end
    disp(savename)
    savename = fullfile(figures_folder,savename);
%     %print2svg(hf, savename, 1, 1);
    
    
    
end % for dc


%% CBF mean +-std, of each donor


% close all

markarr = '^odsv';

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


flag_no72h = 1;
if flag_no72h
    time = 1:3;
else
    time = 1:4;%24 .* (0:3);
end

time_jitter = linspace(-0.1, 0.1, 4);

lfs = 10;   % labels fontsize
tfs = 8;    % ticks fontsize
lw = 1.2;


bsi = find(boxsizes_vector == 64);

% one figure per donor
for dc = 1:numel(donors)
    
    hf = figure;
    hf.Units = 'centimeters';
    hf.Position([3 4]) = [8 6];
    
    ha = axes;
    ha.Box = 'on';
    ha.NextPlot = 'add';
    ha.FontSize = tfs;
    ha.FontName = 'Arial';
    
    clear he
    
    for stc = 1:numel(sampletypes)
        
        
        he(stc) = errorbar(time + time_jitter(stc),...
            arrayfun(@(i) nanmean(MergedData(dc, stc, i).Frequency_Hz{bsi}), time),...
            arrayfun(@(i) nanstd(MergedData(dc, stc, i).Frequency_Hz{bsi}), time));
        he(stc).Color = st_cmap(stc,:);
        he(stc).MarkerFaceColor = 'w';
        he(stc).Marker = markarr(stc);
        he(stc).LineStyle = 'none';
        he(stc).LineWidth = lw;
        
    end %for stc
    
    ha.XTickLabel = {''};
    ha.XTick = time;
    ha.XTickLabel = timepoints(time);
    ha.FontSize = tfs;
    
    ha.YLim = [2 12];
    ha.XLim = [0.5 numel(timepoints(time)) + 0.5];
    
    ha.XLabel.String = 'Time';
    ha.XLabel.FontSize = lfs;
    
    ha.YLabel.String = 'CBF, [Hz]';
    ha.YLabel.FontSize = lfs;
    
    % fix title
    ha.Title.String = strjoin({donors{dc}},' ');
    ha.Title.String = replace(ha.Title.String, {'LumacaFTOR';'Orkambi';'T1a'},{'VX809';'VX809 + VX770';'T\alpha1'});
    ha.Title.FontWeight = 'normal';
    ha.Title.FontSize = lfs;
    ha.Title.String = regexprep(ha.Title.String,'d(?=\d)','s');
    
    hl = legend(he,...
        replace({MergedData(dc,:,1).sampletype_str},...
        {'LumacaFTOR';'Orkambi';'T1a'},{'VX809';'VX809 + VX770';'T\alpha1'}));
    
    
    hl.Box = 'off';
    hl.Orientation = 'vertical';
    hl.FontSize = lfs;
    hl.Units = 'normalized';
    
    ha.Position = [0.1591 0.1546 0.82 0.79];
    hl.Position(1) = ha.Position(1) - 0.02;
    hl.Position(2) = sum(ha.Position([2 4])) - hl.Position(4);
    
    
    savename = strjoin({'CBF_errorbar',...
        [donors{dc}],...
        'alltreatments'},'_');
    if flag_no72h
        savename = [savename,'_no72h'];
    end
    disp(savename)
    savename = fullfile(figures_folder,savename);
%     %print2svg(hf, savename, 1, 1);
    
end %for


%% CBF mean +-std, of each donor divided by the control


% close all

markarr = '^odsv';

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


flag_no72h = 1;
if flag_no72h
    time = 1:3;
else
    time = 1:4;%24 .* (0:3);
end

time_jitter = linspace(-0.1, 0.1, 4);

lfs = 10;   % labels fontsize
tfs = 8;    % ticks fontsize
lw = 1.2;


bsi = find(boxsizes_vector == 64);

% one figure per donor
for dc = 1:numel(donors)
    
    hf = figure;
    hf.Units = 'centimeters';
    hf.Position([3 4]) = [8 6];
    
    ha = axes;
    ha.Box = 'on';
    ha.NextPlot = 'add';
    ha.FontSize = tfs;
    ha.FontName = 'Arial';
    
    clear he
    
    for stc = 1:numel(sampletypes)-1
        
        xplot = time + time_jitter(stc);
        
        yplot = arrayfun(@(i) nanmean(MergedData(dc, stc, i).Frequency_Hz{bsi}), time) ./...
            arrayfun(@(i) nanmean(MergedData(dc, 4, i).Frequency_Hz{bsi}), time);
        
        eplot = hypot(arrayfun(@(i) nanstd(MergedData(dc, stc, i).Frequency_Hz{bsi}) ./ nanmean(MergedData(dc, stc, i).Frequency_Hz{bsi}), time),...
            arrayfun(@(i) nanstd(MergedData(dc, 4, i).Frequency_Hz{bsi}) ./ nanmean(MergedData(dc, 4, i).Frequency_Hz{bsi}), time) ) .* yplot;
        
        
        he(stc) = errorbar(xplot, yplot, eplot);
        he(stc).Color = st_cmap(stc,:);
        he(stc).MarkerFaceColor = 'w';
        he(stc).Marker = markarr(stc);
        he(stc).LineStyle = 'none';
        he(stc).LineWidth = lw;
        
    end %for stc
    
    ha.XTickLabel = {''};
    ha.XTick = time;
    ha.XTickLabel = timepoints(time);
    ha.FontSize = tfs;
    
    ha.YLim = [0 3];
    ha.XLim = [0.5 numel(timepoints(time)) + 0.5];
    
    ha.XLabel.String = 'Time';
    ha.XLabel.FontSize = lfs;
    
    ha.YLabel.String = 'CBF/CBF_{Ctrl}';
    ha.YLabel.FontSize = lfs;
    
    % fix title
    ha.Title.String = strjoin({donors{dc}},' ');
    ha.Title.String = replace(ha.Title.String, {'LumacaFTOR';'Orkambi';'T1a'},{'VX809';'VX809 + VX770';'T\alpha1'});
    ha.Title.FontWeight = 'normal';
    ha.Title.FontSize = lfs;
    ha.Title.String = regexprep(ha.Title.String,'d(?=\d)','s');
    
    hl = legend(he,...
        replace({MergedData(dc,1:end-1,1).sampletype_str},...
        {'LumacaFTOR';'Orkambi';'T1a'},{'VX809';'VX809 + VX770';'T\alpha1'}));
    
    
    hl.Box = 'off';
    hl.Orientation = 'vertical';
    hl.FontSize = lfs;
    hl.Units = 'normalized';
    
    ha.Position = [0.1591 0.1546 0.82 0.79];
    hl.Position(1) = ha.Position(1) - 0.02;
    hl.Position(2) = sum(ha.Position([2 4])) - hl.Position(4);
    
    
    savename = strjoin({'CBF_errorbar',...
        [donors{dc}],...
        'alltreatments_CtrlNormalised'},'_');
    if flag_no72h
        savename = [savename,'_no72h'];
    end
    disp(savename)
    savename = fullfile(figures_folder,savename);
%     %print2svg(hf, savename, 1, 1);
    
end %for



%% CBF(t) - CBF_ctrl(t) mean +-std, of each donor 


% close all

markarr = '^odsv';

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


flag_no72h = 1;
if flag_no72h
    time = 1:3;
else
    time = 1:4;%24 .* (0:3);
end

time_jitter = linspace(-0.1, 0.1, 4);

lfs = 10;   % labels fontsize
tfs = 8;    % ticks fontsize
lw = 1.2;


bsi = find(boxsizes_vector == 64);

% one figure per donor
for dc = 1:numel(donors)
    
    hf = figure;
    hf.Units = 'centimeters';
    hf.Position([3 4]) = [8 6];
    
    ha = axes;
    ha.Box = 'on';
    ha.NextPlot = 'add';
    ha.FontSize = tfs;
    ha.FontName = 'Arial';
    
    clear he
    
    for stc = 1:numel(sampletypes)-1
        
        xplot = time + time_jitter(stc);
        
        yplot = arrayfun(@(i) nanmean(MergedData(dc, stc, i).Frequency_Hz{bsi}), time) - ...
            arrayfun(@(i)nanmean(MergedData(dc, 4, i).Frequency_Hz{bsi}), time) ;
        
        eplot = hypot(arrayfun(@(i) nanstd(MergedData(dc, stc, i).Frequency_Hz{bsi}), time) , ...
            arrayfun(@(i)nanstd(MergedData(dc, 4, i).Frequency_Hz{bsi}), time) );
        
%         eplot(1) = 0; %number minus itself is 0
        
        
        he(stc) = errorbar(xplot, yplot, eplot);
        he(stc).Color = st_cmap(stc,:);
        he(stc).MarkerFaceColor = 'w';
        he(stc).Marker = markarr(stc);
        he(stc).LineStyle = 'none';
        he(stc).LineWidth = lw;
        
    end %for stc
    
    ha.XTickLabel = {''};
    ha.XTick = time;
    ha.XTickLabel = timepoints(time);
    ha.FontSize = tfs;
    
    ha.YLim = [-4 6];
    ha.XLim = [0.5 numel(timepoints(time)) + 0.5];
    
    ha.XLabel.String = 'Time';
    ha.XLabel.FontSize = lfs;
    
    ha.YLabel.String = '\DeltaCBF(t), [Hz]';
    ha.YLabel.FontSize = lfs;
    
    % fix title
    ha.Title.String = strjoin({donors{dc}},' ');
    ha.Title.String = replace(ha.Title.String, {'LumacaFTOR';'Orkambi';'T1a'},{'VX809';'VX809 + VX770';'T\alpha1'});
    ha.Title.FontWeight = 'normal';
    ha.Title.FontSize = lfs;
    ha.Title.String = regexprep(ha.Title.String,'d(?=\d)','s');
    
    hl = legend(he,...
        replace({MergedData(dc,:,1).sampletype_str},...
        {'LumacaFTOR';'Orkambi';'T1a'},{'VX809';'VX809 + VX770';'T\alpha1'}));
    
    
    hl.Box = 'off';
    hl.Orientation = 'vertical';
    hl.FontSize = lfs;
    hl.Units = 'normalized';
    
    ha.Position = [0.1591 0.1546 0.82 0.79];
    hl.Position(1) = ha.Position(1) - 0.02;
    hl.Position(2) = sum(ha.Position([2 4])) - hl.Position(4);
    
    
    savename = strjoin({'DeltaCBF_errorbar',...
        [donors{dc}],...
        'alltreatments'},'_');
    if flag_no72h
        savename = [savename,'_no72h'];
    end
    disp(savename)
    savename = fullfile(figures_folder,savename);
%     %print2svg(hf, savename, 1, 1);
    
end %for

%% CBF errorbar taking mean of (median of each donor)
% STILL TODO




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
        hl.String = regexprep(hl.String,'d(?=\d)','s');
        
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
        %print2svg(hf, savename, 1, 1);
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

shiftedmuMat = reshape(shiftedmuMat,size(MergedData));
confint_shiftedmuMat = reshape(confint_shiftedmuMat,size(MergedData));
timeMat = reshape(timeMat,size(MergedData));

plottimeMat = 1+timeMat./24 + repmat(linspace(-0.1, 0.1, 4),3,1,4);


st_cmap = vertcat(darkB, darkG, darkO, darkGrey) ;

flag_no72h = 1;

if flag_no72h
    ind_plot = 1:3;
else
    ind_plot = 1:4;
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
for dc = 1:3
    
    
    % plot data
    he(dc,:) = plot(squeeze(plottimeMat(dc,:,ind_plot))',...
        10.^squeeze(shiftedmuMat(dc,:,ind_plot))');
      
    
    % appearances
    for stc = 1:4
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

% legend
hl = legend(he(1:3:end),...
    replace({MergedData(1,:,1).sampletype_str}, {'LumacaFTOR';'Orkambi';'T1a'},{'Lumacaftor';'Lumacaftor + Ivacaftor';'T\alpha1'}));


hl.Box = 'off';
hl.Orientation = 'vertical';
hl.FontSize = lfs;
hl.Units = 'normalized';

% axes properties
ha.XTick = 1:size(MergedData,2);
ha.XLim = minmax(ind_plot) + 0.5.*[-1 1];
ha.XTickLabel = {MergedData(1,1,ind_plot).timepoint_str} ;
% ha.XAxis.FontSize = lfs;
ha.XLabel.String = 'Time';
ha.XLabel.FontSize = lfs;
ha.YLabel.String = '\lambda^2, [\mum^2]';
ha.YLabel.FontSize = lfs;
ha.YLim = [1e2 10^7.1];
setsemilogy

% positions
ha.Position = [0.1591 0.1546 0.82 0.79];
hl.Position(1) = ha.Position(1) - 0.02;
hl.Position(2) = sum(ha.Position([2 4])) - hl.Position(4);
ha.YLabel.Units = 'normalized';
ha.YLabel.Position(1) = -0.13;
ha.YLabel.VerticalAlignment = 'baseline';

savename = fullfile('left_lambda2(t)');
if flag_no72h
    savename = [savename,'_no72h'];
end
savename = fullfile(figures_folder,savename);
% %print2svg(hf, savename, 1, 1);

disp('log10lambda(0h), one row per drug')
disp(shiftedmuMat(:,:,1)')
disp('err log10lambda(0h), one row per drug')
disp(2*confint_shiftedmuMat(:,:,1)')
disp('log10lambda(48h), one row per drug')
disp(shiftedmuMat(:,:,3)')
disp('err log10lambda(48h), one row per drug')
disp(2*confint_shiftedmuMat(:,:,3)')


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
flag_trendline = 0;
flag_no72h = 1;

if flag_no72h
    ind_plot = 1:3;
else
    ind_plot = 1:4;
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


%% normalise each donavg timepoint with the donavg control at the same time

ctrl_norm_shiftedmu_donavg = shiftedmu_donavg(1:3,:) - shiftedmu_donavg(4,:);
err_ctrl_shiftednorm_mu_donavg = sqrt( shiftedmu_err(1:3,:).^2 + shiftedmu_err(end,:).^2  );

% %%close all

flag_trendline = 1;
flag_no72h = 1;

if flag_no72h
    ind_plot = 1:3;
else
    ind_plot = 1:4;
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
% exponentiate
y = 10.^ctrl_norm_shiftedmu_donavg(:,ind_plot)';
yneg = y - 10.^(ctrl_norm_shiftedmu_donavg(:,ind_plot)' - err_ctrl_shiftednorm_mu_donavg(:,ind_plot)');
ypos = 10.^(ctrl_norm_shiftedmu_donavg(:,ind_plot)' + err_ctrl_shiftednorm_mu_donavg(:,ind_plot)') - y;
x = (repmat(ind_plot,3,1) + 0.12.*[-1;0;1])';

he = errorbar(x,y,yneg,ypos,'o');

% legend
hl = legend(he, replace({MergedData(1,1:end-1,1).sampletype_str},...
    {'LumacaFTOR';'Orkambi';'T1a'},{'Lumacaftor';'Lumacaftor + Ivacaftor';'T\alpha1'}));

for stc = 1:numel(sampletypes)-1
    
    % appearances
    he(stc).Marker = 'o';
    he(stc).MarkerFaceColor = 'w';
    he(stc).LineWidth = 1.2;
    he(stc).CapSize = 5;
    % set(he_m,'LineStyle','none');
    he(stc).Color = st_cmap(stc,:);
    
    % create trendline
    if flag_trendline
        fitline = fit(ind_plot', ctrl_norm_shiftedmu_donavg(stc,ind_plot)','poly1',...
            fitoptions('weights',1./err_ctrl_shiftednorm_mu_donavg(stc,ind_plot)));
        
        hpf(stc) = plot(x(ind_plot,stc), 10.^fitline(ind_plot) );
        hpf(stc).Marker = 'none';
        hpf(stc).LineWidth = 1.2;
        hpf(stc).Color = st_cmap(stc,:);
    end %if
    
end %for stc

% axes properties
ha.FontSize = 8;
ha.XTick = ind_plot;
ha.XLim = minmax(ha.XTick) + 0.5.*[-1 1];
ha.XTickLabel = {MergedData(1,1,ind_plot).timepoint_str} ;

ha.XLabel.String = 'Time';
ha.XLabel.FontSize = lfs;
ha.YLabel.String = '<\lambda^2(t)>_d / <\lambda^2_{Ctrl}(t)>_d';
ha.YLabel.FontSize = lfs;



ha.YLim = [1e-2 1e1];
setsemilogy

hl.Box = 'off';
hl.FontSize = lfs;
hl.Units = 'normalized';

% positions
ha.Position = [0.1591 0.1546 0.82 0.79];
hl.Position(1) = ha.Position(1) - 0.01;
hl.Position(2) = ha.Position(2) + 0.01;
ha.YLabel.Units = 'normalized';
ha.YLabel.Position(1) = -0.13;
ha.YLabel.VerticalAlignment = 'baseline';


savename = 'left_lambda2(t)_average_CtrlNormalised';

if flag_trendline
    savename = [savename,'_trendline'];
end

if flag_no72h
    savename = [savename,'_no72h'];
end
savename = fullfile(figures_folder,savename);
%print2svg(hf, savename, 1, 1);

%% now plot log10(lambda2(48h)) - log10(lambda2(0h)) for each donor/treatment

%%close all

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
st_cmap = vertcat(darkB, darkG, darkO, darkGrey) ;

x_stagger = [-0.09 0 0.09];

lfs = 10;   % labels fontsize
tfs = 8;    % ticks fontsize
lw = 1.2;

hf = figure;
hf.Units = 'centimeters';
hf.Position([3 4]) = [8 6];

ha = axes;
ha.NextPlot = 'add';
ha.Box = 'on';
ha.FontName = 'Arial';

clear he;

matforpvalues = nan(3,4);

% sampletype on x axes, 
for stc = 1:numel(sampletypes)
    for dc = 1:numel(donors)
        
        % no need to apply the left shift as we are then subtracting one
        % from the other
        mu_48h = MergedData(dc,stc,3).Damping_Hz_fit_out2.mu;
        mu_48h_ci = diff(par_confint(MergedData(dc,stc,3).Damping_Hz_fit_out2,'mu',0.68))./2;
        mu_00h = MergedData(dc,stc,1).Damping_Hz_fit_out2.mu;
        mu_00h_ci = diff(par_confint(MergedData(dc,stc,1).Damping_Hz_fit_out2,'mu',0.68))./2;
        
        delta_mu = mu_48h - mu_00h;
        err_delta_mu = sqrt(mu_48h_ci.^2 + mu_00h_ci.^2);
        
        delta_lambda2 = 10.^delta_mu;
        ue_delta_lambda2 = 10.^(delta_mu + err_delta_mu) - 10.^(delta_mu);
        le_delta_lambda2 = 10.^(delta_mu) - 10.^(delta_mu - err_delta_mu);
                
%         he(stc,dc) = errorbar(stc + x_stagger(dc), delta_mu, err_delta_mu);
        he(stc,dc) = errorbar(stc + x_stagger(dc), delta_lambda2,...
            le_delta_lambda2, ue_delta_lambda2);

        
        he(stc,dc).Marker = markarr(dc);
        he(stc,dc).MarkerFaceColor = 'w';
        he(stc,dc).Color = st_cmap(stc,:);
        he(stc,dc).CapSize = 5;
        he(stc,dc).LineWidth = lw;
        he(stc,dc).LineStyle = 'none';
        
        matforpvalues(dc,stc) = delta_mu;
        errmatforpvalues(dc,stc) = 2*err_delta_mu; %twice because err_delta_mu was halved for plotting
        
    end %for dc
end %for

setsemilogy
ha.XLim = [0.3 4.3];


hl = legend(he(end,:),donors);
hl.Box = 'off';
hl.FontSize = lfs;
hl.Units = 'normalized';
hl.String = regexprep(hl.String,'d(?=\d)','s');

% axes properties
ha.FontSize = 8;
ha.XTick = 1:numel(sampletypes);
ha.XLim = minmax(ha.XTick) + 0.5.*[-1 1];
ha.XTickLabel = sampletypes ;
ha.XTickLabel = replace(ha.XTickLabel,...
     {'LumacaFTOR';'Orkambi';'T1a'},{'Lumacaftor';'Ivacaftor Lumacaftor';'T\alpha1'});
ha.XAxis.FontSize = lfs;
% ha.XLabel.String = ' ';
% ha.XLabel.FontSize = lfs;
ha.YLabel.String = '\lambda^2_{48h} / \lambda^2_{00h}';
ha.YLabel.FontSize = lfs;
ht = fix_xticklabels(ha);
set(ht,'FontSize',lfs,...
    'FontName','Arial',...
    'Color',darkGrey);


% positions
ha.Position = [0.1591 0.1546 0.82 0.79];
hl.Position(1) = ha.Position(1) - 0.01;
hl.Position(2) = ha.Position(2) + 0.01;
ha.YLabel.Units = 'normalized';
ha.YLabel.Position(1) = -0.13;
ha.YLabel.VerticalAlignment = 'baseline';

savename = 'lambda2_ratio_48h';
savename = fullfile(figures_folder,savename);
%print2svg(hf, savename, 1, 1);

disp('using paired t-test')
for stc = 1:numel(sampletypes)-1
    disp(sampletypes{stc})
    disp('P value \lambda2(48h) / \lambda2(0h) vs control');
    [~,p] = ttest(matforpvalues(:,stc), matforpvalues(:,end),'tail','left');
    disp(p)
end

disp('using unpaired t-test')
for stc = 1:numel(sampletypes)-1
    disp(sampletypes{stc})
    disp('P value \lambda2(48h) / \lambda2(0h) vs control');
    [~,p] = ttest2(matforpvalues(:,stc), matforpvalues(:,end),'tail','left');
    disp(p)
end

matforpvalues
errmatforpvalues



%% now plot log10(lambda2(48h)) - log10(lambda2(0h)) - [ log10(lambda2(48h)) - log10(lambda2(0h))]ctrl
% for each donor/treatment

%%close all

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
st_cmap = vertcat(darkB, darkG, darkO, darkGrey) ;

x_stagger = [-0.09 0 0.09];

lfs = 10;   % labels fontsize
tfs = 8;    % ticks fontsize
lw = 1.2;

hf = figure;
hf.Units = 'centimeters';
hf.Position([3 4]) = [8 6];

ha = axes;
ha.NextPlot = 'add';
ha.Box = 'on';
ha.FontName = 'Arial';

clear he;

matforpvalues = nan(3,4);
errmatforpvalues = nan(3,4);

% sampletype on x axes, 
for stc = 1:numel(sampletypes)-1
    for dc = 1:numel(donors)
        
        % no need to apply the left shift as we are then subtracting one
        % from the other
        mu_48h = MergedData(dc,stc,3).Damping_Hz_fit_out2.mu;
        mu_48h_ci = diff(par_confint(MergedData(dc,stc,3).Damping_Hz_fit_out2,'mu',0.68))./2;
        mu_00h = MergedData(dc,stc,1).Damping_Hz_fit_out2.mu;
        mu_00h_ci = diff(par_confint(MergedData(dc,stc,1).Damping_Hz_fit_out2,'mu',0.68))./2;
        
        delta_mu = mu_48h - mu_00h;
        err_delta_mu = sqrt(mu_48h_ci.^2 + mu_00h_ci.^2);
        
        % now same for control (yeah I shouldn;t recalculate it every time
        % but whatever)
        mu_48h_ctrl = MergedData(dc,end,3).Damping_Hz_fit_out2.mu;
        mu_48h_ci_ctrl = diff(par_confint(MergedData(dc,end,3).Damping_Hz_fit_out2,'mu',0.68))./2;
        mu_00h_ctrl = MergedData(dc,end,1).Damping_Hz_fit_out2.mu;
        mu_00h_ci_ctrl = diff(par_confint(MergedData(dc,end,1).Damping_Hz_fit_out2,'mu',0.68))./2;
        
        delta_mu_ctrl = mu_48h_ctrl - mu_00h_ctrl;
        err_delta_mu_ctrl = sqrt(mu_48h_ci_ctrl.^2 + mu_00h_ci_ctrl.^2);
        
        delta2_mu = delta_mu - delta_mu_ctrl;
        err_delta2_mu = sqrt(err_delta_mu.^2 + err_delta_mu_ctrl.^2);
        
        delta2_lambda2 = 10.^delta2_mu;
        10.^delta_mu ./ 10.^delta_mu_ctrl;
        ue_delta2_lambda2 = 10.^(delta2_mu + err_delta2_mu) - 10.^(delta2_mu);
        le_delta2_lambda2 = 10.^(delta2_mu) - 10.^(delta2_mu - err_delta2_mu);
                
%         he(stc,dc) = errorbar(stc + x_stagger(dc), delta_mu, err_delta_mu);
        he(stc,dc) = errorbar(stc + x_stagger(dc), delta2_lambda2,...
            le_delta2_lambda2, ue_delta2_lambda2);

        
        he(stc,dc).Marker = markarr(dc);
        he(stc,dc).MarkerFaceColor = 'w';
        he(stc,dc).Color = st_cmap(stc,:);
        he(stc,dc).CapSize = 5;
        he(stc,dc).LineWidth = lw;
        he(stc,dc).LineStyle = 'none';
        
        matforpvalues(dc,stc) = delta2_mu;
        errmatforpvalues(dc,stc) = 2*err_delta2_mu; %twice because err_delta_mu was halved for plotting
        
    end %for dc
end %for

setsemilogy
ha.XLim = [0.3 4.3];


hl = legend(he(end,:),donors);
hl.Box = 'off';
hl.FontSize = lfs;
hl.Units = 'normalized';
hl.String = regexprep(hl.String,'d(?=\d)','s');

% axes properties
ha.FontSize = 8;
ha.XTick = 1:numel(sampletypes)-1;
ha.XLim = minmax(ha.XTick) + 0.5.*[-1 1];
ha.XTickLabel = sampletypes ;
ha.XTickLabel = replace(ha.XTickLabel,...
    {'LumacaFTOR';'Orkambi';'T1a'},{'VX809';'[VX809+VX770]';'T\alpha1'});
%      {'LumacaFTOR';'Orkambi';'T1a'},{'Lumacaftor';'Ivacaftor Lumacaftor';'T\alpha1'});
ha.XAxis.FontSize = lfs;
% ha.XLabel.String = ' ';
% ha.XLabel.FontSize = lfs;
ha.YLabel.String = '\lambda^2_{48h} / \lambda^2_{00h} / (\lambda^2_{48h} / \lambda^2_{00h})_{ctrl}';
ha.YLabel.FontSize = lfs;
% ht = fix_xticklabels(ha);
set(ht,'FontSize',lfs,...
    'FontName','Arial',...
    'Color',darkGrey);


% positions
ha.Position = [0.1591 0.1546 0.82 0.79];
hl.Position(1) = ha.Position(1) - 0.01;
hl.Position(2) = ha.Position(2) + 0.01;
ha.YLabel.Units = 'normalized';
ha.YLabel.Position(1) = -0.13;
ha.YLabel.VerticalAlignment = 'baseline';

plot(xlim,ones(2,1),'--','color',[0.5 0.5 0.5],'linewidth',0.4)

savename = 'lambda2_ratio_48h_normalised_by_control';
savename = fullfile(figures_folder,savename);
print2svg(hf, savename, 1, 1);

disp('using paired t-test')
for stc = 1:numel(sampletypes)-1
    disp(sampletypes{stc})
    disp('P value \lambda2(48h) / \lambda2(0h) vs control');
    [~,p] = ttest(matforpvalues(:,stc), 0,'tail','left');
    disp(p)
end



matforpvalues
errmatforpvalues



%% now plot lamda2/lambda2(0h) for each donor/treatment

%%close all


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
st_cmap = vertcat(darkB, darkG, darkO, darkGrey) ;

% x_stagger = [-0.09 0 0.09];
x_stagger = linspace(-0.1, 0.1,4);

flag_no72h = 1;

if flag_no72h
    ind_plot = 1:3;
else
    ind_plot = 1:4;
end

lfs = 10;   % labels fontsize
tfs = 8;    % ticks fontsize
lw = 1.2;

hf = figure;
hf.Units = 'centimeters';
hf.Position([3 4]) = [8 6];

ha = axes;
ha.NextPlot = 'add';
ha.Box = 'on';
ha.FontName = 'Arial';

clear he;

% sampletype on x axes, 
for stc = 1:numel(sampletypes)
    for dc = 1:numel(donors)
        
        % no need to apply the left shift as we are then subtracting one
        % from the other
        mu_t = arrayfun(@(i)MergedData(dc,stc,i).Damping_Hz_fit_out2.mu,ind_plot);
        mu_t_ci = arrayfun(@(i)diff(par_confint(MergedData(dc,stc,i).Damping_Hz_fit_out2,'mu',0.68))./2,ind_plot);
        
        nmu_t = mu_t - mu_t(1);
        nmu_t_ci = sqrt( mu_t_ci.^2 + mu_t_ci(1).^2 );
        
                
        norm_lambda2 = 10.^nmu_t;
        ue_norm_lambda2 = 10.^(nmu_t + nmu_t_ci) - 10.^(nmu_t);
        le_norm_lambda2 = 10.^(nmu_t) - 10.^(nmu_t - nmu_t_ci);
                
        he(stc,dc) = plot(ind_plot + x_stagger(stc), norm_lambda2);

        
        he(stc,dc).Marker = markarr(dc);
        he(stc,dc).MarkerSize = 6;
        he(stc,dc).MarkerFaceColor = 'w';
        he(stc,dc).Color = st_cmap(stc,:);
%         he(stc,dc).CapSize = 5;
        he(stc,dc).LineWidth = lw;
        he(stc,dc).LineStyle = 'none';
        
    end %for dc
end %for

setsemilogy


hl = legend(he(:,1),sampletypes);
hl.String = replace(hl.String,{'LumacaFTOR';'Orkambi';'T1a'},{'Lumacaftor';'Ivacaftor + Lumacaftor';'T\alpha1'});
hl.Box = 'off';
hl.FontSize = lfs;
hl.Units = 'normalized';

% axes properties
ha.FontSize = 8;
ha.XTick = ind_plot;
ha.XLim = minmax(ha.XTick) + 0.5.*[-1 1];
ha.XTickLabel = timepoints(ind_plot);

ha.YLim = [1e-2 1e5];

ha.XLabel.String = 'Time';
ha.XLabel.FontSize = lfs;
ha.YLabel.String = '\lambda^2(t) / \lambda^2(0)';
ha.YLabel.FontSize = lfs;



% positions
ha.Position = [0.1591 0.1546 0.82 0.79];
hl.Position(1) = ha.Position(1) - 0.02;
hl.Position(2) = sum(ha.Position([2 4])) - hl.Position(4);
ha.YLabel.Units = 'normalized';
ha.YLabel.Position(1) = -0.13;
ha.YLabel.VerticalAlignment = 'baseline';

savename = 'lambda2_ratio(t)';
if flag_no72h
    savename = [savename,'_no72h'];
end
savename = fullfile(figures_folder,savename);
%print2svg(hf, savename, 1, 1);


%% same but averaging across donors

%%close all


% divide each lambda by its 00h point (so subtract mu_00h from mu(t)
norm_mu_t = shiftedmuMat - shiftedmuMat(:,:,1); 
norm_mu_t_ci = sqrt( confint_shiftedmuMat.^2 + confint_shiftedmuMat(:,:,1).^2 ); 
% set the errorbar of time 0 at 0 by definition (dividing a number by
% itself)
norm_mu_t_ci(:,:,1) = eps;
% (this is basically what was plotted in the previous figure)

% now average across donors
[norm_mu_t_donavg, err_norm_mu_t_donavg] = weightednanmean(norm_mu_t,norm_mu_t_ci,1);
norm_mu_t_donavg = squeeze(norm_mu_t_donavg);
err_norm_mu_t_donavg = squeeze(err_norm_mu_t_donavg);
% now time is 2nd dimension, sampletype is 1st

flag_trendline = 0;
flag_no72h = 1;

if flag_no72h
    ind_plot = 1:3;
else
    ind_plot = 1:4;
end

hf = figure;
hf.Units = 'centimeters';
hf.Position([3 4]) = [8 6];

ha = axes;
ha.NextPlot = 'add';
ha.Box = 'on';
ha.FontName = 'Arial';

clear he;

y = 10.^norm_mu_t_donavg';
yneg = y - 10.^(norm_mu_t_donavg - err_norm_mu_t_donavg)';
ypos = 10.^(norm_mu_t_donavg + err_norm_mu_t_donavg)' - y;
x = (repmat(1:4,4,1) + linspace(-0.05,0.05,4)')';

he = errorbar(x(ind_plot,:),y(ind_plot,:),yneg(ind_plot,:),ypos(ind_plot,:));


for stc = 1:numel(sampletypes)
    
    % appearances
    he(stc).Marker = 'o';
    he(stc).MarkerFaceColor = 'w';
    he(stc).LineWidth = 1.2;
    he(stc).CapSize = 5;
    he(stc).LineStyle = 'none';
    he(stc).Color = st_cmap(stc,:);
    
    % create trendline
    if flag_trendline
        fitline = fit(ind_plot', norm_mu_t_donavg(stc,ind_plot)','poly1',...
            fitoptions('weights',1./err_norm_mu_t_donavg(stc,ind_plot)));
        
        hpf(stc) = plot(x(ind_plot,stc), 10.^fitline(ind_plot) );
        hpf(stc).Marker = 'none';
        hpf(stc).LineWidth = 1.2;
        hpf(stc).Color = st_cmap(stc,:);
    end %if
%     
end %for stc


setsemilogy


hl = legend(he,sampletypes);
hl.String = replace(hl.String,{'LumacaFTOR';'Orkambi';'T1a'},{'Lumacaftor';'Ivacaftor + Lumacaftor';'T\alpha1'});

hl.Box = 'off';
hl.FontSize = lfs;
hl.Units = 'normalized';

% axes properties
ha.FontSize = 8;
ha.XTick = ind_plot;
ha.XLim = minmax(ha.XTick) + 0.5.*[-1 1];
ha.XTickLabel = timepoints(ind_plot);

ha.YLim = [1e-2 1e1];

ha.XLabel.String = 'Time';
ha.XLabel.FontSize = lfs;
ha.YLabel.String = '<\lambda^2(t) / \lambda^2(0)>_d';
ha.YLabel.FontSize = lfs;

% positions
ha.Position = [0.1591 0.1546 0.82 0.79];
hl.Position(1) = ha.Position(1) - 0.01;
hl.Position(2) = ha.Position(2) + 0.01;
ha.YLabel.Units = 'normalized';
ha.YLabel.Position(1) = -0.13;
ha.YLabel.VerticalAlignment = 'baseline';


savename = 'lambda2_ratio(t)_averaged';
if flag_trendline
    savename = [savename,'_trendline'];
end

if flag_no72h
    savename = [savename,'_no72h'];
end
savename = fullfile(figures_folder,savename);
%print2svg(hf, savename, 1, 1);

%% <lambda(t)>_d / <lambda(0)>_d Don't know if we want this

%%close all

% average by donor
[mu_t_donavg, err_mu_t_donavg] = weightednanmean(shiftedmuMat,confint_shiftedmuMat,1);
mu_t_donavg = squeeze(mu_t_donavg);
err_mu_t_donavg = squeeze(err_mu_t_donavg);
% now time is 2nd dimension, sampletype is 1st



% divide each <lambda(t)>d by its 00h point (so subtract mu_00h from mu(t)
norm_mu_t_donavg = mu_t_donavg - mu_t_donavg(:,1); 
err_norm_mu_t_donavg = sqrt( err_mu_t_donavg.^2 + err_mu_t_donavg(:,1).^2 ); 

% set the errorbar of time 0 at 0 by definition (dividing a number by
% itself)
err_norm_mu_t_donavg(:,1) = eps;


flag_trendline = 1;
flag_no72h = 1;

if flag_no72h
    ind_plot = 1:3;
else
    ind_plot = 1:4;
end

hf = figure;
hf.Units = 'centimeters';
hf.Position([3 4]) = [8 6];

ha = axes;
ha.NextPlot = 'add';
ha.Box = 'on';
ha.FontName = 'Arial';

clear he;

y = 10.^norm_mu_t_donavg';
yneg = y - 10.^(norm_mu_t_donavg - err_norm_mu_t_donavg)';
ypos = 10.^(norm_mu_t_donavg + err_norm_mu_t_donavg)' - y;
x = (repmat(1:4,4,1) + linspace(-0.05,0.05,4)')';

he = errorbar(x(ind_plot,:),y(ind_plot,:),yneg(ind_plot,:),ypos(ind_plot,:));


for stc = 1:numel(sampletypes)
    
    % appearances
    he(stc).Marker = 'o';
    he(stc).MarkerFaceColor = 'w';
    he(stc).LineWidth = 1.2;
    he(stc).CapSize = 5;
    he(stc).LineStyle = 'none';
    he(stc).Color = st_cmap(stc,:);
    
    % create trendline
    if flag_trendline
        fitline = fit(ind_plot', norm_mu_t_donavg(stc,ind_plot)','poly1',...
            fitoptions('weights',1./err_norm_mu_t_donavg(stc,ind_plot)));
        
        hpf(stc) = plot(x(ind_plot,stc), 10.^fitline(ind_plot) );
        hpf(stc).Marker = 'none';
        hpf(stc).LineWidth = 1.2;
        hpf(stc).Color = st_cmap(stc,:);
    end %if
%     
end %for stc


setsemilogy


hl = legend(he,sampletypes);
hl.String = replace(hl.String,{'LumacaFTOR';'Orkambi';'T1a'},{'Lumacaftor';'Ivacaftor + Lumacaftor';'T\alpha1'});

hl.Box = 'off';
hl.FontSize = lfs;
hl.Units = 'normalized';

% axes properties
ha.FontSize = 8;
ha.XTick = ind_plot;
ha.XLim = minmax(ha.XTick) + 0.5.*[-1 1];
ha.XTickLabel = timepoints(ind_plot);

ha.YLim = [1e-2 1e1];

ha.XLabel.String = 'Time';
ha.XLabel.FontSize = lfs;
ha.YLabel.String = '<\lambda^2(t)>_d / <\lambda^2(0)>_d';
ha.YLabel.FontSize = lfs;

% positions
ha.Position = [0.1591 0.1546 0.82 0.79];
hl.Position(1) = ha.Position(1) - 0.01;
hl.Position(2) = ha.Position(2) + 0.01;
ha.YLabel.Units = 'normalized';
ha.YLabel.Position(1) = -0.13;
ha.YLabel.VerticalAlignment = 'baseline';


% savename = 'lambda2_ratio(t)_averaged';
% if flag_trendline
%     savename = [savename,'_trendline'];
% end
% 
% if flag_no72h
%     savename = [savename,'_no72h'];
% end
% __%print2svg(hf,savename,1,1);


%% normalise each timeseries by its time 0, then normalise by the control as well

%%close all


% divide each lambda by its 00h point (so subtract mu_00h from mu(t)
norm_mu_t = shiftedmuMat - shiftedmuMat(:,:,1); 
norm_mu_t_ci = sqrt( confint_shiftedmuMat.^2 + confint_shiftedmuMat(:,:,1).^2 ); 
% set the errorbar of time 0 at 0 by definition (dividing a number by
% itself)
norm_mu_t_ci(:,:,1) = eps;

% now divide each donor by the corespondent control
% well subtract as it's a log thing
norm2_mu_t = norm_mu_t(:,1:numel(sampletypes)-1,:) - norm_mu_t(:,numel(sampletypes),:);% name bc normalised twice lol
norm2_mu_t_ci = sqrt( norm_mu_t_ci(:,1:numel(sampletypes)-1,:).^2 + norm_mu_t_ci(:,numel(sampletypes),:).^2 ); 
% "error" by propagating conf ints

flag_no72h = 1;

if flag_no72h
    ind_plot = 1:3;
else
    ind_plot = 1:4;
end

hf = figure;
hf.Units = 'centimeters';
hf.Position([3 4]) = [8 6];

ha = axes;
ha.NextPlot = 'add';
ha.Box = 'on';
ha.FontName = 'Arial';

clear he;

% sampletype on x axes, 
for stc = 1:numel(sampletypes)-1
    for dc = 1:numel(donors)
        
     
        y = 10.^norm2_mu_t(dc,stc,ind_plot);
        yneg = y - 10.^(norm2_mu_t(dc,stc,ind_plot) - norm2_mu_t_ci(dc,stc,ind_plot));
        ypos = 10.^(norm2_mu_t(dc,stc,ind_plot) + norm2_mu_t_ci(dc,stc,ind_plot)) - y;
        y = squeeze(y);
        yneg = squeeze(yneg);
        ypos = squeeze(ypos);
        
        he(stc,dc) = plot(ind_plot + x_stagger(stc), y);

        he(stc,dc).Marker = markarr(dc);
        he(stc,dc).MarkerSize = 6;
        he(stc,dc).MarkerFaceColor = 'w';
        he(stc,dc).Color = st_cmap(stc,:);
%         he(stc,dc).CapSize = 5;
        he(stc,dc).LineWidth = lw;
        he(stc,dc).LineStyle = 'none';
        
    end %for dc
end %for

setsemilogy


hl = legend(he(:,1),sampletypes(1:end-1));
hl.String = replace(hl.String,{'LumacaFTOR';'Orkambi';'T1a'},{'Lumacaftor';'Ivacaftor + Lumacaftor';'T\alpha1'});

hl.Box = 'off';
hl.FontSize = lfs;
hl.Units = 'normalized';

% axes properties
ha.FontSize = 8;
ha.XTick = ind_plot;
ha.XLim = minmax(ha.XTick) + 0.5.*[-1 1];
ha.XTickLabel = timepoints(ind_plot);

ha.YLim = [1e-4 1e1];

ha.XLabel.String = 'Time';
ha.XLabel.FontSize = lfs;
ha.YLabel.String = '(\lambda^2 / \lambda^2(0)) / (\lambda^2_{Ctrl} / \lambda^2_{Ctrl}(0))';
ha.YLabel.FontSize = lfs;

% positions
ha.Position = [0.1591 0.1546 0.82 0.79];
hl.Position(1) = ha.Position(1) - 0.01;
hl.Position(2) = ha.Position(2) + 0.01;
ha.YLabel.Units = 'normalized';
ha.YLabel.Position(1) = -0.13;
ha.YLabel.VerticalAlignment = 'baseline';


savename = 'lambda2_ratio(t)_CtrlNormalised';

if flag_no72h
    savename = [savename,'_no72h'];
end
savename = fullfile(figures_folder,savename);
%print2svg(hf, savename, 1, 1);


%% now same quantity but average across donors

[norm2_mu_t_donavg, err_norm2_mu_t_donavg] = weightednanmean(norm2_mu_t,norm2_mu_t_ci,1);
norm2_mu_t_donavg = squeeze(norm2_mu_t_donavg);
err_norm2_mu_t_donavg = squeeze(err_norm2_mu_t_donavg);
% now time is 2nd dimension, sampletype is 1st, and first errorbar is ~eps
% as dividing 1 by 1 anyway

flag_trendline = 0;
flag_no72h = 1;

if flag_no72h
    ind_plot = 1:3;
else
    ind_plot = 1:4;
end

x = ind_plot' + [-0.05 0 0.05];
y = 10.^norm2_mu_t_donavg(:,ind_plot)';
yneg = y - 10.^(norm2_mu_t_donavg(:,ind_plot) - err_norm2_mu_t_donavg(:,ind_plot))';
ypos = 10.^(norm2_mu_t_donavg(:,ind_plot) + err_norm2_mu_t_donavg(:,ind_plot))' - y;


%%close all

hf = figure;
hf.Units = 'centimeters';
hf.Position([3 4]) = [8 6];

ha = axes;
ha.NextPlot = 'add';
ha.Box = 'on';
ha.FontName = 'Arial';

he = errorbar(x, y, yneg, ypos);

for stc = 1:numel(sampletypes) - 1
    
    % appearances
    he(stc).Marker = 'o';
    he(stc).MarkerFaceColor = 'w';
    he(stc).LineWidth = 1.2;
    he(stc).CapSize = 5;
    he(stc).LineStyle = 'none';
    he(stc).Color = st_cmap(stc,:);
    
    % create trendline
    if flag_trendline
        fitline = fit(ind_plot', norm2_mu_t_donavg(stc,ind_plot)','poly1',...
            fitoptions('weights',1./err_norm2_mu_t_donavg(stc,ind_plot)));
        
        hpf(stc) = plot(x(ind_plot,stc), 10.^fitline(ind_plot) );
        hpf(stc).Marker = 'none';
        hpf(stc).LineWidth = 1.2;
        hpf(stc).Color = st_cmap(stc,:);
    end %if
%     
end %for stc
 


setsemilogy


hl = legend(he,sampletypes(1:end-1));
hl.String = replace(hl.String,{'LumacaFTOR';'Orkambi';'T1a'},{'Lumacaftor';'Ivacaftor + Lumacaftor';'T\alpha1'});

hl.Box = 'off';
hl.FontSize = lfs;
hl.Units = 'normalized';

% axes properties
ha.FontSize = 8;
ha.XTick = ind_plot;
ha.XLim = minmax(ha.XTick) + 0.5.*[-1 1];
ha.XTickLabel = timepoints(ind_plot);

ha.YLim = [1e-1 1e1];

ha.XLabel.String = 'Time';
ha.XLabel.FontSize = lfs;
ha.YLabel.String = '<(\lambda^2 / \lambda^2(0)) / (\lambda^2_{Ctrl} / \lambda^2_{Ctrl}(0))>_d';
ha.YLabel.FontSize = lfs;

% positions
ha.Position = [0.1591 0.1546 0.82 0.79];
hl.Position(1) = ha.Position(1) - 0.01;
hl.Position(2) = ha.Position(2) + 0.01;
ha.YLabel.Units = 'normalized';
ha.YLabel.Position(1) = -0.13;
ha.YLabel.VerticalAlignment = 'baseline';


savename = 'lambda2_ratio(t)_CtrlNormalised_average';
if flag_trendline
    savename = [savename,'_trendline'];
end

if flag_no72h
    savename = [savename,'_no72h'];
end
savename = fullfile(figures_folder,savename);
%print2svg(hf, savename, 1, 1);



