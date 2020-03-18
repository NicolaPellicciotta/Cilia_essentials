function [hf, ha, he, hpf, hpv, hpp] = plot_sigmoid_errorbar_compare(MergedData, colormapname, ha, flag_lambda2error, flag_errorbar, whichshoulder)
%plot_sigmoid_errorbar_compare takes a MergedData structure or an array of them (as output from
%merge_SampleType_data) and plots sigmoid with the fitting line
%and errorbar

%% input check

% colormap
if nargin < 2 || isempty(colormapname) || ~ischar(colormapname) || ~exist(colormapname,'file')
    if numel(MergedData) <= 3
        colormapname = 'lines';
    else 
        colormapname = 'parula';
    end %if
end

if nargin < 3 || isempty(ha)
    % no axes in input, create figure and axes
    hf = figure;
    ha = axes;
else
    hf = ha.Parent;
end %if

if nargin < 4 || isempty(flag_lambda2error)
    flag_lambda2error = false;
end

if nargin < 5 || isempty(flag_errorbar)
    flag_errorbar = true;
end

if nargin < 6 || isempty(whichshoulder) || ~ischar(whichshoulder) || ~any(strcmpi(whichshoulder,{'left','right'}))
    whichshoulder = 'right';
end %if

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

ha.Box = 'on';
ha.NextPlot = 'add';

% axes limits
setsemilogx
ha.YLim = [0 max( arrayfun(@(i) max(vertcat(MergedData(i).Damping_Hz{:})), 1:numel(MergedData)) )];
ha.XLim = [1e0 1e5];
ha.XTick = logspace(0,5,6);

% plots
for i = 1:numel(MergedData)

	% errorbars
    if flag_errorbar
	he(i) = errorbar(MergedData(i).window_area_um2, MergedData(i).med_Damping_Hz,...
        MergedData(i).ler_Damping_Hz, MergedData(i).uer_Damping_Hz);
    else
        he(i) = plot(MergedData(i).window_area_um2, MergedData(i).med_Damping_Hz);
    end
	he(i).LineWidth = 1.2;
	he(i).LineStyle = 'none';
	he(i).Marker = 'o';
    he(i).MarkerSize = 6;
	he(i).MarkerFaceColor = 'w';
	he(i).MarkerEdgeColor = cmap(i,:);
	he(i).Color = cmap(i,:);
    
	% fit
% 	hpf(i) = plot(xx, MergedData(i).Damping_Hz_fit_out(xx));
    hpf(i) = plot(xx, MergedData(i).Damping_Hz_fit_out2(log10(xx)));    
	hpf(i).LineWidth = 1.2;
	hpf(i).Color = he(i).MarkerEdgeColor;
	
	% vertical line
    switch whichshoulder
        case 'right'
%             hpv(i) = plot(MergedData(i).Damping_Hz_fit_out.b*[1 1]*exp(2), [0 20*ha.YLim(2)]);
            hpv(i) = plot(10.^(MergedData(i).Damping_Hz_fit_out2.mu)*[1 1].*exp(2), [0 20*ha.YLim(2)]);
        case 'left'
%             hpv(i) = plot(MergedData(i).Damping_Hz_fit_out.b*[1 1]/exp(2), [0 20*ha.YLim(2)]);
            hpv(i) = plot(10^(MergedData(i).Damping_Hz_fit_out2.mu)*[1 1]./exp(2), [0 20*ha.YLim(2)]);
        otherwise
    end
	hpv(i).Color = hpf(i).Color;
	hpv(i).LineStyle = hpf(i).LineStyle;
	hpv(i).LineWidth = hpf(i).LineWidth;
    

    if flag_lambda2error
%         dummy = par_confint(MergedData(i).Damping_Hz_fit_out,'b');
        dummy = 10.^(par_confint(MergedData(i).Damping_Hz_fit_out2,'mu',0.68));
        
        switch whichshoulder
            case 'right'
                hpp(i) = patch('XData',dummy([1 2 2 1]).*exp(2),'YData', ha.YLim([1 1 2 2]));
            case 'left'
                hpp(i) = patch('XData',dummy([1 2 2 1])./exp(2),'YData', ha.YLim([1 1 2 2]));
            otherwise
        end
        hpp(i).FaceColor = hpf(i).Color;
        hpp(i).FaceAlpha = 0.2;
        hpp(i).EdgeColor = 'none';
        

        
%         hplv(i) = plot(dummy(1,2)*[1 1], [0 20*ha.YLim(2)]);
%         hplv(i).Color = 
%         hplv(i).LineStyle = '--';
%         hplv(i).LineWidth = hpf(i).LineWidth;
%         
%         hprv(i) = plot(dummy(2,2)*[1 1], [0 20*ha.YLim(2)]);
%         hprv(i).Color = hpf(i).Color;
%         hprv(i).LineStyle = '--';
%         hprv(i).LineWidth = hpf(i).LineWidth;
    end
    
    % prepare legend entry
    dummy_ca = arrayfun(@(ff)strjoin(cellstr(MergedData(i).(ff{:}))),...
        fnames_leg,'UniformOutput',0);
    leg{i} = strjoin(dummy_ca,', ');
end


% legend
hleg = legend(ha, he, leg);
hleg.Box = 'off';
hleg.Location = 'NorthWest';
hleg.FontSize = 10;
hleg.Interpreter = 'none';


% labels and titles
ha.XLabel.String = 'DDM Window Area, [\mum^2]';
ha.YLabel.String = '1/\tau_c, [s^{-1}]';
ha.XLabel.FontSize = 16;
ha.YLabel.FontSize = 16;


