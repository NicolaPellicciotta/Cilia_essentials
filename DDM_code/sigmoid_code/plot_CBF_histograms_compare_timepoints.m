function [ hf ] = plot_CBF_histograms_compare_timepoints( MergedData, BoxSizeIndex, binedges, colormapname, flag_numbers )
%plot_CBF_histograms_compare_timepoints Takes the structure array MergedData
%and plots the CBF in a one-axis histogram. if there's more than one
%timepoint, these will go on the same axis in different colors
%   MergedData can either be a scalar structure from merge_SampleType_data
%   or it can be an array of such structure. I'm assuming that the
%   sampletype does not change, and that the only thigs that changes
%   between entries of the structure array is the timepoints


%% input check

% first check that the conditions for MergedData have been met

% check that sampletype_str does not change
if any(~arrayfun(@(i)strcmp(MergedData(1).sampletype_str, MergedData(i).sampletype_str), 1:numel(MergedData)) )
    error('Can only plot same sample type');
end

if any(~cell2mat(arrayfun(@(i)strcmp(MergedData(1).inserts_str, MergedData(i).inserts_str),...
        1:numel(MergedData), 'UniformOutput', 0)) )
    error('Can only plot same inserts');
end

if any(~cell2mat(arrayfun(@(i)strcmp(MergedData(1).donors_str, MergedData(i).donors_str),...
        1:numel(MergedData), 'UniformOutput', 0)) )
    error('Can only plot same donors');
end

% now the rest of the input check

% flag numbers
if nargin < 5 || isempty(flag_numbers)
    flag_numbers = false;
end

% colormap
if nargin < 4 || isempty(colormapname) || ~ischar(colormapname) || ~exist(colormapname,'file')
    if numel(MergedData) <= 3
        colormapname = 'lines';
    else 
        colormapname = 'parula';
    end %if
end


% binedges
if nargin < 3 || isempty(binedges) || ~isvector(binedges)
    binedges = 0:1:30;
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

cmap = eval([colormapname,'(',num2str(numel(MergedData)),')']);

%% plot

% prepare figure
hf = figure;

% prepare axes
ha = axes;
% ha.Position = [0.15 0.15 0.8 0.75];
ha.XLim = minmax(binedges);
ha.Box = 'on';
ha.NextPlot = 'add';

% labels
ha.XLabel.String = 'CBF, [Hz]';
ha.XLabel.FontSize = 16;
ha.YLabel.String = 'Counts, normalized';
ha.YLabel.FontSize = 16;
ha.Title.String = [sprintf('%s    ',MergedData(1).sampletype_str),sprintf('%s    ',MergedData(1).donors_str{:})];
ha.Title.Interpreter = 'none';
ha.Title.FontSize = 16;

% actual plot
for tpc = 1:numel(MergedData)
    
    % histogram with edge
    hhs(tpc) = histogram(ha, MergedData(tpc).Frequency_Hz{bsi}, binedges);
    hhs(tpc).Normalization = 'Probability';
    hhs(tpc).DisplayStyle = 'Stairs';
    hhs(tpc).LineWidth = 1.5;
    hhs(tpc).EdgeColor = cmap(tpc,:);
    
    % edgeless histogram
    hhb(tpc) = histogram(ha, MergedData(tpc).Frequency_Hz{bsi}, binedges);
    hhb(tpc).Normalization = 'Probability';
    hhb(tpc).EdgeColor = 'none';
    hhb(tpc).FaceColor = cmap(tpc,:);
    hhb(tpc).FaceAlpha = 0.2;

    % legend
    if flag_numbers
        leg{tpc} = [MergedData(tpc).timepoint_str, ', (',...
            num2str(nanmean(MergedData(tpc).Frequency_Hz{bsi}),3),...
            ' ± ',num2str(nanstd(MergedData(tpc).Frequency_Hz{bsi}),3), ')Hz'];
    else
        leg{tpc} = MergedData(tpc).timepoint_str;
    end %if
    
end %for


% legend
hleg = legend(ha, hhb, leg);
hleg.Box = 'off';
hleg.FontSize = 10;
hleg.Interpreter = 'none';



end

