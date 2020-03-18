function [hs, ha3, ha, ht1, hrb1] = Iqtau_reader( Iqtau, fits )
%Iqtau_reader simple GUI to scroll through the modes of an Iqtau matrix


%% input check

if nargin < 2
    fits = [];
end %if

%% set parameters

[max_mode, max_lag] = size(Iqtau);
frame_size = max_mode * 2;


%% set up figure

screensize = get(0,'ScreenSize');

hf = figure;
hf.Position = [(screensize(3) - hf.Position(3)*2)/2 hf.Position(2) hf.Position(3)*2 hf.Position(4)];


%% Initial figure with mode

ha3 = axes;
ha3.Units = 'normalized';
ha3.Position = [0.13/2 0.1227 0.74/2 0.7912];
ha3.Box = 'on';
ha3.NextPlot = 'replacechildren';
% ha3.XScale = 'Log';

ha3.XLabel.String = 'Lag time, [frames]';
ha3.XLabel.FontSize = 16;
ha3.XLabel.Units = 'normalized';
ha3xlpos = ha3.XLabel.Position;
ha3.YLabel.String = sprintf('I(q=%.2f,\tau), [a.u.]',pi*1/max_mode);
ha3.YLabel.FontSize = 16;
ha3.YLabel.Units = 'normalized';
ha3.YLabel.Position = [-0.07 0.5 0];

%% Initial figure with Iqtau

ha = axes('Units','normalized');
ha.Position = [0.5 + 0.13/2 0.1227 0.74/2 0.7912];

hi = imagesc(Iqtau);

ha.YLabel.String = 'Fourier Mode';
ha.YLabel.FontSize = 16;
ha.XLabel.String = 'Lag time, [frames]';
ha.XLabel.FontSize = 16;
ha.Title.String = 'Image Structure Function I(q,\tau)';
ha.Title.FontSize = 18;
ha.NextPlot = 'add';

% colorbar
hc = colorbar;
hc.Units = 'normalized';
hc.Location = 'South';
hc.Color = 'w';
hc.Label.String = 'I(q,\tau), [a.u.]';
hc.Label.FontSize = 14;

% copy axes to put alternative units
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

% fix positions
ha2.Position = ha.Position;

% better background color, and prepare for printing
hf.Color = 'w';
hf.PaperPositionMode = 'auto';
hf.InvertHardcopy = 'off';

%% set slider

hs = uicontrol('Style','Slider');
hs.Min = 1;
hs.Max = max_mode;
hs.Value = max_mode;
hs.Units = 'Normalized';
hs.Position = [sum(ha3.Position([1 3]))+0.13/4 ha3.Position(2)+0.1 0.13/4 ha3.Position(4)-0.1];
hs.BackgroundColor = 'w';
hs.SliderStep = [1 10]./max_mode;

hl = addlistener(hs,'Value','PostSet',@update_figure);

%% set text (fixed)

ht1 = uicontrol('Style','Text');
ht1.String = {'Fourier'; 'Mode'};
ht1.FontSize = 12;
ht1.Units = 'normalized';
ht1.Position = [hs.Position(1)-0.02 ha3.Position(2)-0.01 0.13/4+0.04 0.1];
ht1.BackgroundColor = 'w';

%% set text (editable)

ht2 = uicontrol('Style','Text');
ht2.FontSize = 12;
ht2.Units = 'normalized';
ht2.Position = [hs.Position(1)-0.02 ha3.Position(2)-0.04 0.13/4+0.04 0.04];
ht2.BackgroundColor = 'w';

%% set linlin - semilogx radio button

hbg = uibuttongroup(hf);
hbg.Title = 'XScale';
hbg.Units = 'normalized';
hbg.Position = [ ha3.Position(1) + ha3.Position(3)/2 - 0.06  sum(ha3.Position([2 4])) 0.12 0.08];
hbg.BackgroundColor = 'w';

hrb1 = uicontrol('Style','radiobutton');
hrb1.Units = 'normalized';
hrb1.Position = [hbg.Position(1) + 0.1 * hbg.Position(3), hbg.Position(2) + 0.1*hbg.Position(4), [0.3,0.4].*hbg.Position([3 4])];
hrb1.String = 'lin';
hrb1.BackgroundColor = 'w';
hrb1.Value = true;
hrb1.Callback = @hrb1_callback;

hrb2 = uicontrol('Style','radiobutton');
hrb2.Units = 'normalized';
hrb2.Position = [hbg.Position(1) + 0.6 * hbg.Position(3), hbg.Position(2) + 0.1*hbg.Position(4), [0.3,0.4].*hbg.Position([3 4])];
hrb2.String = 'log';
hrb2.BackgroundColor = 'w';
hrb2.Value = false;
hrb2.Callback = @hrb2_callback;

%% plot first mode
update_figure;

%% update figure

    function [] = update_figure(~,~)
        
        mode_to_show = max_mode - round(hs.Value) + 1;
        
        % update plot
        axes(ha3)
        ha3.YLabel.String = sprintf('I(q=%.2f px^{-1},\\tau), [a.u.]',pi*mode_to_show/max_mode);
        
        % and plot fits if they are there
        if ~isempty(fits)
            
            hp = plot(ha3,1:max_lag, Iqtau(mode_to_show,1:max_lag)./fits(mode_to_show).norm_fact,...
            '.','MarkerSize',8);
            ha3.NextPlot = 'add';
            xplot = 0:.1:max_lag;
            plot(xplot, fits(mode_to_show).fit_out(xplot),...
                'r-','LineWidth',1.2);
            ha3.NextPlot = 'replacechildren';
            
        else % no fits, no normalisation either
            
            hp = plot(ha3,1:max_lag, Iqtau(mode_to_show,1:max_lag),...
            '.','MarkerSize',8);
        
        end %if
        
        % fix position of xlabel
        ha3.XLabel.Position = ha3xlpos;
        
        % update red line on imagesc
        
        axes(ha)
        delete(findobj(ha,'Type','Line'));
        hrl = plot(ha.XLim, mode_to_show*[1 1], 'r', 'LineWidth', 1.2);
        
        % update text
        ht2.String = num2str(mode_to_show);
        
    end %function


    function [] = hrb1_callback(~,~)
        
        ha3.XScale = 'lin';
        hrb2.Value = false;
        update_figure;
        
    end %function
    
    function [] = hrb2_callback(~,~)
        
        ha3.XScale = 'log';
        hrb1.Value = false;
        update_figure;
        
    end %function
end