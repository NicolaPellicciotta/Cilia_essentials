function [] = plot_single_sigmoid_CBFhistograms(MergedData)
%plot_single_sigmoid_CBFhistograms takes a MergedData structure (as output from
%merge_SampleType_data) and plots a check sigmoid with all the CBF
%histograms (to check that it does not depend too much on the box size)

% prepare figure
hf = figure;

% prepare axes
NBoxsizes = numel(MergedData.window_area_um2);
for bsc = 1:NBoxsizes
    ha(bsc) = axes;
    ha(bsc).Position = [0.1+(bsc-1)*0.8/NBoxsizes, 0.15, 0.8/NBoxsizes, 0.75];
    ha(bsc).YLim = [0 max(vertcat(MergedData.Frequency_Hz{:}))];
%     ha(bsc).XLim = [0 1];
    ha(bsc).Box = 'on';
    ha(bsc).NextPlot = 'add';
    if bsc > 1
        ha(bsc).YTickLabels = [];
    end %if
end %for



% histograms and median
for bsc = 1:NBoxsizes
    hh(bsc) = histogram(ha(bsc),MergedData.Frequency_Hz{bsc},...
        'Normalization','Probability',...
        'Orientation','Horizontal',...
        'DisplayStyle','Stairs',...
        'LineWidth',2);
    hhd(bsc) = histogram(ha(bsc),MergedData.Frequency_Hz{bsc},...
        'Normalization','Probability',...
        'Orientation','Horizontal',...
        'EdgeColor','none',...
        'FaceColor', ha(bsc).ColorOrder(1,:),...
        'LineWidth',2);
    hpm(bsc) = plot(ha(bsc), ha(bsc).XLim, median(MergedData.Frequency_Hz{bsc})*[1 1]);
    hpm(bsc).LineWidth = 1.2;
    hpm(bsc).Color = ha(bsc).ColorOrder(2,:);
end %for

 
% labels and titles
for bsc = 1:NBoxsizes
    ha(bsc).Title.String = sprintf('%4.2f \\mum^2',MergedData.window_area_um2(bsc));
    ha(bsc).Title.FontWeight = 'normal';
end %for
ha(1).YLabel.String = 'CBF, [Hz]';
ha(1).YLabel.FontSize = 16;

% dummy axes for labels
hda = axes;
hda.Position(1) = min(horzcat(arrayfun(@(i)ha(i).Position(1), 1:NBoxsizes)));
hda.Position(2) = ha(1).Position(2);
hda.Position(3) = sum(horzcat(arrayfun(@(i)ha(i).Position(3), 1:NBoxsizes)));
hda.Position(4) = ha(1).Position(4);
hda.Visible = 'off';
hda.XLabel.Visible = 'on';
hda.XLabel.String = 'Counts, normalised';
hda.XLabel.FontSize = 16;




