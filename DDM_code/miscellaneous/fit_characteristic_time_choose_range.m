function [hf, min_mode_fit, max_mode_fit] = fit_characteristic_time_choose_range(modes, tau_c, err_tau_c)
% make the user choose the mode range to fit

%% initial parameters, modified via GUI
min_mode_fit_temp = modes(1);
max_mode_fit_temp = modes(end);

%% figure parameters

hf = figure('CloseRequestFcn',@my_closereq);

ha = axes;
ha.Box = 'on';
hold on;


%% plot parameters

he = errorbar(modes, tau_c, err_tau_c);
he.Marker = '.';
he.MarkerSize = 8;
he.Color = 'b';
he.LineStyle = 'none';
he.LineWidth = 1.2;

%% labels

hxl = xlabel('Fourier Modes');
hxl.FontSize = 16;

hyl = ylabel('\tau_c, [frames]');
hyl.FontSize = 16;

ht = title({'Move the sliders to select the range'; 'you want to fit, then close this figure'});
ht.FontSize = 16;
ht.FontWeight = 'bold';

%% set loglog now because errorbar disrupts it

ha.XScale = 'log';
ha.YScale = 'log';

drawnow

%% sliders

% move a bit the figure first
ha.Position(2) = 0.3;
ha.Position(4) = 0.55;

hs1 = uicontrol('Style','Slider');
hs1.Min = min_mode_fit_temp;
hs1.Max = max_mode_fit_temp;
hs1.Value = min_mode_fit_temp;
hs1.Units = 'Normalized';
hs1.Position = [ha.Position(1) 0.1 ha.Position(3) 0.04];
hs1.BackgroundColor = 'w';
hs1.SliderStep = [1 10]./(max_mode_fit_temp-min_mode_fit_temp +1);

hl1 = addlistener(hs1,'Value','PostSet',@update_figure);

hs2 = uicontrol('Style','Slider');
hs2.Min = min_mode_fit_temp;
hs2.Max = max_mode_fit_temp;
hs2.Value = max_mode_fit_temp;
hs2.Units = 'Normalized';
hs2.Position = [ha.Position(1) 0.05 ha.Position(3) 0.04];
hs2.BackgroundColor = 'w';
hs2.SliderStep = [1 10]./(max_mode_fit_temp-min_mode_fit_temp +1);

hl2 = addlistener(hs2,'Value','PostSet',@update_figure);

%% update

update_figure;

    function [] = update_figure(~,~)
        
        min_mode_fit_temp = min(round(hs1.Value), round(hs2.Value));
        max_mode_fit_temp = max(round(hs1.Value), round(hs2.Value));
        
        hprevplot = findobj(ha,'Type','Errorbar','Color','g');
        if ~isempty(hprevplot)
            delete(hprevplot)
        end
        
        he1 = errorbar(modes(min_mode_fit_temp:max_mode_fit_temp),...
            tau_c(min_mode_fit_temp:max_mode_fit_temp),...
            err_tau_c(min_mode_fit_temp:max_mode_fit_temp));
        he1.Marker = '.';
        he1.MarkerSize = 8;
        he1.Color = 'g';
        he1.LineStyle = 'none';
        he1.LineWidth = 1.2;
        
        % set again scale
        ha.XScale = 'log';
        ha.YScale = 'log';
        
        % push out value
        assignin('caller','min_mode_fit',min_mode_fit_temp);
        assignin('caller','max_mode_fit',max_mode_fit_temp);
        
        drawnow
        
    end %function

    function [] = my_closereq(~,~)
        
        min_mode_fit = min(round(hs1.Value), round(hs2.Value));
        max_mode_fit = max(round(hs1.Value), round(hs2.Value));
        
        
        delete(gcf)
    end
end
