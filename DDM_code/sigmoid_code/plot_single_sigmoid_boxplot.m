function [] = plot_single_sigmoid_boxplot(MergedData)
%plot_single_sigmoid_boxplot takes a MergedData structure (as output from
%merge_SampleType_data) and plots a single sigmoid using boxplots

% prepare figure
hf = figure;

% prepare axes
ha = axes;
box on;
hold on;

% axes limits
% setsemilogx
% ha.YLim = [0 max(vertcat(MergedData.Damping_Hz{:}))];
% ha.XLim = [1e1 1e5];
% ha.XTick = logspace(1,5,5);

% prepare data for boxplots
NBoxsizes = numel(MergedData.window_area_um2);
XX = vertcat(MergedData.Damping_Hz{:});
GG = vertcat(cell2mat( arrayfun( @(i)repmat(MergedData.window_area_um2(i), numel(MergedData.Damping_Hz{i}),1),...
    (1:NBoxsizes)', 'UniformOutput',false ) ));

% boxplot
hb = boxplot(XX, GG, 'Position', log(MergedData.window_area_um2));


 
% labels and titles
ha.XLabel.String = 'DDM Window Area, [\mum^2]';
ha.YLabel.String = '1/\tau_c, [s^{-1}]';
ha.XLabel.FontSize = 16;
ha.YLabel.FontSize = 16;


