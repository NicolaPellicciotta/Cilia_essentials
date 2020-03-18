function [ hf, ha ] = plot_CBF_histograms_compare_sampletypes( MergedData, BoxSizeIndex, binedges, colormapname, flag_numbers )
%plot_CBF_histograms_compare_sampletypes Takes the structure array MergedData
%and plots the CBF in multiple axes histograms. if there's more than one
%timepoint, these will go on different axes.
%   I'm assuming that what changes is sampletype and timepoint, and that
%   MergedData is a matrix in which the first dimension is the sample type,
%   the second dimension is timepoints


%% input check

% first check that the conditions for MergedData have been met

% check that sampletype_str does not change along each row
if any(~arrayfun(@(i,j)strcmp(MergedData(i,1).sampletype_str, MergedData(i,j).sampletype_str),...
    repmat((1:size(MergedData,1))',1,size(MergedData,2)),...
    repmat(1:size(MergedData,2), size(MergedData,1), 1) ))
    error('Each row has to have the same sampletype');
end

% check that timepoint_str does not change along each column
if any(~arrayfun(@(i,j)strcmp(MergedData(1,j).timepoint_str, MergedData(i,j).timepoint_str),...
    repmat((1:size(MergedData,1))',1,size(MergedData,2)),...
    repmat(1:size(MergedData,2), size(MergedData,1), 1) ))
    error('Each column has to have the same timepoint_str');
end

% check insert is the same along a row
if any(~cell2mat(arrayfun(@(i,j)strcmp(MergedData(i,1).inserts_str, MergedData(i,j).inserts_str),...
        repmat((1:size(MergedData,1))',1,size(MergedData,2)),...
        repmat(1:size(MergedData,2), size(MergedData,1), 1),...
        'UniformOutput', 0)) )
    error('Can only plot same inserts');
end

% check donor is the same along a row
if any(~cell2mat(arrayfun(@(i,j)strcmp(MergedData(i,1).donors_str, MergedData(i,j).donors_str),...
        repmat((1:size(MergedData,1))',1,size(MergedData,2)),...
        repmat(1:size(MergedData,2), size(MergedData,1), 1),...
        'UniformOutput', 0)) )
    error('Can only plot same donors');
end

% now the rest of the input check

% flag numbers
if nargin < 5 || isempty(flag_numbers)
    flag_numbers = false;
end

% colormap
if nargin < 4 || isempty(colormapname) || ~ischar(colormapname) || ~exist(colormapname,'file')
    if size(MergedData,1) <= 3
        colormapname = 'lines';
    else
        colormapname = 'parula';
    end %if
end


% binedges
if nargin < 3 || isempty(binedges) || ~isvector(binedges)
    binedges = 0:1:45;
end %if


% boxsize index
if nargin < 2 || isempty(BoxSizeIndex) || BoxSizeIndex > numel(MergedData(1).window_area_um2)
    bsi = min(3, numel(MergedData(1).window_area_um2) );
else
    bsi = BoxSizeIndex;
end %if



%% automatic process fix parameters

binedges = sort(binedges);
binedges(binedges < 0 | binedges > 60) = [];

cmap = eval([colormapname,'(',num2str(size(MergedData,1)),')']);

%% plot

% how many axes? as many as timepoints
nax = size(MergedData,2);
xloff = 0.1;
xroff = 0.05;
yboff = 0.15;
ytoff = 0.1;
graph_width = (1-xloff - xroff)/nax;
graph_height = 1 - yboff - ytoff;

% prepare figure
hf = figure;
if nax > 1
    centred_resized_figure(hf,[1 nax/2]);
end %if

for tpc = 1:size(MergedData,2)
    
    % prepare axes
    ha(tpc) = axes;
    ha(tpc).Position = [xloff + (tpc-1)*graph_width, yboff, graph_width, graph_height ];
    ha(tpc).YLim = minmax(binedges);
    ha(tpc).Box = 'on';
    ha(tpc).NextPlot = 'add';
    
    % labels
    ha(tpc).XLabel.String = 'Counts, normalized';
    ha(tpc).XLabel.FontSize = 16;
    if tpc == 1
        ha(tpc).YLabel.String = 'CBF, [Hz]';
        ha(tpc).YLabel.FontSize = 16;
    else
        ha(tpc).YTickLabels = [];
    end
%     ha(tpc).Title.String = [sprintf('%s    ',MergedData(1).sampletype_str),sprintf('%s    ',MergedData(1).donors_str{:})];
    ha(tpc).Title.String = [sprintf('%s    ',MergedData(1,tpc).timepoint_str),sprintf('%s    ',MergedData(1, tpc).donors_str{:})];
    ha(tpc).Title.Interpreter = 'none';
    ha(tpc).Title.FontSize = 16;
    
    % actual plot
    for stc = 1:size(MergedData,1)
        % histogram with edge
        hhs(stc, tpc) = histogram(ha(tpc), MergedData(stc, tpc).Frequency_Hz{bsi}, binedges);
        hhs(stc, tpc).Normalization = 'Probability';
        hhs(stc, tpc).DisplayStyle = 'Stairs';
        hhs(stc, tpc).LineWidth = 1.5;
        hhs(stc, tpc).EdgeColor = cmap(stc,:);
        hhs(stc, tpc).Orientation = 'Horizontal';
        
        % edgeless histogram
        hhb(stc, tpc) = histogram(ha(tpc), MergedData(stc, tpc).Frequency_Hz{bsi}, binedges);
        hhb(stc, tpc).Normalization = 'Probability';
        hhb(stc, tpc).EdgeColor = 'none';
        hhb(stc, tpc).FaceColor = cmap(stc,:);
        hhb(stc, tpc).FaceAlpha = 0.2;
        hhb(stc, tpc).Orientation = 'Horizontal';
        
        % legend
        if flag_numbers
            leg{stc, tpc} = [MergedData(stc, tpc).sampletype_str, ', (',...
                num2str(nanmean(MergedData(stc, tpc).Frequency_Hz{bsi}),3),...
                ' � ',num2str(nanstd(MergedData(stc, tpc).Frequency_Hz{bsi}),3), ')Hz'];
        else
            leg{stc, tpc} = MergedData(stc, tpc).sampletype_str;
        end %if
        
    end %for
    
    
    % legend
    hleg(tpc) = legend(ha(tpc), hhb(:,tpc), leg{:,tpc});
    hleg(tpc).Box = 'off';
    hleg(tpc).FontSize = 10;
    hleg(tpc).Interpreter = 'none';
    
    
    
end


% set univocuous x lims
set(ha, 'XLim', [0 max(reshape(cell2mat(get(ha,'XLim')),[],1))])

% for tpc = 1:size(MergedData, 2)
%     ha(tpc).XTickLabel{end} = '';
% end


