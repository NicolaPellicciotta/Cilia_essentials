function [hf, ha] = plot_sigmoid_all_errorbar_compare_sampletypes(MergedData, colormapname)
%plot_sigmoid_errorbar_compare takes a MergedData structure or an array of them (as output from
%merge_SampleType_data) and plots sigmoid with the fitting line
%and errorbar
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

if any(~cell2mat(arrayfun(@(i)strcmp(MergedData(1).inserts_str, MergedData(i).inserts_str),...
        1:numel(MergedData), 'UniformOutput', 0)) )
    error('Can only plot same inserts');
end

if any(~cell2mat(arrayfun(@(i)strcmp(MergedData(1).donors_str, MergedData(i).donors_str),...
        1:numel(MergedData), 'UniformOutput', 0)) )
    error('Can only plot same donors');
end

% colormap
if nargin < 2 || isempty(colormapname) || ~ischar(colormapname) || ~exist(colormapname,'file')
    if numel(MergedData) <= 3
        colormapname = 'lines';
    else
        colormapname = 'parula';
    end %if
end


%% parameters

% prepare x for fit
xx = logspace(0,5,1e3);

% colormap
cmap = eval([colormapname,'(',num2str(numel(MergedData)),')']);


%% find which strings differ between mergeddata entries

fnames = fieldnames(MergedData);
fnames = fnames(contains(fnames,'_str'));

idx_changing_fnames = ~arrayfun(@(i)isequal(MergedData(:).(i{:})), fnames);
fnames_leg = fnames(idx_changing_fnames);

%% plot

% prepare figure
naxc = size(MergedData,2);
naxr = 2;
xloff = 0.1;
xroff = 0.05;
yboff = 0.15;
ytoff = 0.1;
ysep = 0.2;
graph_width = (1-xloff - xroff)/naxc;
graph_height = (1 - yboff - ytoff- ysep)/naxr;

% prepare figure
hf = figure;
if naxc > 1
    centred_resized_figure(hf,[2 naxc/1.5]);
end %if


for tpc = 1:size(MergedData,2)
    
    % ----- first row, not normalised
    
    % prepare axes
    ha(1, tpc) = axes;
    ha(1, tpc).Position = [xloff + (tpc-1)*graph_width, yboff + ysep + graph_height,...
        graph_width, graph_height ];
    
    % and plot on them
    %[~, ha(1, tpc)] = plot_sigmoid_errorbar_compare(MergedData(:,tpc), 'lines', ha(1, tpc), 1);
    [~, ha(1, tpc)] = plot_sigmoid_errorbar_compare(MergedData(:,tpc), 'lines', ha(1, tpc), 1,[],'left');
    % make title
    ha(1, tpc).Title.String = [sprintf('%s    ',MergedData(1,tpc).timepoint_str),sprintf('%s    ',MergedData(1, tpc).donors_str{:})];
    ha(1, tpc).Title.Interpreter = 'none';
    ha(1, tpc).Title.FontSize = 16;
    
    
    
    % ----- second row, normalised

    % prepare axes
    ha(2, tpc) = axes;
    ha(2, tpc).Position = [xloff + (tpc-1)*graph_width, yboff,...
        graph_width, graph_height ];
    
    % and plot on them
    [~, ha(2, tpc)] = plot_sigmoid_normalised_errorbar_compare(MergedData(:,tpc), 'lines', ha(2, tpc), 1);
    
    % make title
    ha(2, tpc).Title.String = [sprintf('%s    ',MergedData(1,tpc).timepoint_str),sprintf('%s    ',MergedData(1, tpc).donors_str{:})];
    ha(2, tpc).Title.Interpreter = 'none';
    ha(2, tpc).Title.FontSize = 16;
    
    
end %for tpc


% now make all axes on a row use the same limits
for i = 1:2
    set(ha(i,:), 'YLim', minmax(cell2mat(get(ha(i,:),'YLim'))) );
    
    % and delete yticklabels and ylabels
    for tpc = 2:size(MergedData,2)
        ha(i, tpc).YTickLabels = [];
        ha(i, tpc).YLabel.String = '';
    end
    
    % remove last ticklabel
    for tpc = 1:naxc-1
        ha(i,tpc).XTickLabels{end} = '';
    end

end %for i












