function [hf, ha, ha2, hc] = Iqtau_show(Iqtau)

%% get parameters

[max_mode, max_lag] = size(Iqtau);
frame_size = max_mode * 2;


%% set handles and plot

hf = figure;
ha = axes('Units','normalized');
hi = surf(Iqtau);


%% set labels

ha.YLabel.String = 'Fourier Mode';
ha.YLabel.FontSize = 16;

ha.XLabel.String = 'Lag time, [frames]';
ha.XLabel.FontSize = 16;

ha.Title.String = 'Image Structure Function I(q,\tau)';
ha.Title.FontSize = 18;


%% set ColorBar

hc = colorbar;

hc.Units = 'normalized';

hc.Location = 'South';

hc.Color = 'w';

hc.Label.String = 'I(q,\tau), [a.u.]';
hc.Label.FontSize = 14;



%% copy axes to put alternative units

ha2 = axes;
ha2.Units = 'normalized';

ha2.Color = 'none';

ha2.YAxisLocation = 'right';
ha2.YDir = 'reverse';

ha2.XLim = ha.XLim;
ha2.XTick = [];

ha2.YLim = ha.YLim;
ha2.YTick = ha.XTick;
ha2.YTickLabel = sprintf('%.1f\n',pi*get(ha2,'YTick')/max_mode);

ha2.YLabel.String = 'Spatial Frequency q, [px^{-1}]';
ha2.YLabel.FontSize = 16;


%% fix positions

ha.Position(3) = 1-2*ha.Position(1);
ha2.Position = ha.Position;


%% better background color, and prepare for printing

hf.Color = 'w';
hf.PaperPositionMode = 'auto';
hf.InvertHardcopy = 'off';


end