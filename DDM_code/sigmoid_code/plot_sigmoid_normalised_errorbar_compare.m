function [hf, ha] = plot_sigmoid_normalised_errorbar_compare(MergedData, colormapname, ha, flag_error)
%plot_sigmoid_normalised_errorbar_compare takes a MergedData structure or an array of them (as output from
%merge_SampleType_data) and plots normalised sigmoid with the fitting line
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


if nargin < 4 || isempty(flag_error)
    flag_error = false;
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

ha.Box = 'on';
ha.NextPlot = 'add';

% axes limits
setsemilogx
ha.YLim = [0 3];
ha.XLim = [1e0 1e5];
ha.XTick = logspace(0,5,6);

% plots
for i = 1:numel(MergedData)

	% errorbars
	he(i) = errorbar(MergedData(i).window_area_um2, MergedData(i).med_Damping_1ocycles,...
    MergedData(i).ler_Damping_1ocycles, MergedData(i).uer_Damping_1ocycles);
	he(i).LineWidth = 1.2;
	he(i).LineStyle = 'none';
	he(i).Marker = 'o';
	he(i).MarkerFaceColor = 'w';
	he(i).MarkerEdgeColor = cmap(i,:);
	he(i).Color = cmap(i,:);
    
	% fit
	hpf(i) = plot(xx, MergedData(i).Damping_1ocycles_fit_out(xx));
	hpf(i).LineWidth = 1.2;
	hpf(i).Color = 0.6 .* he(i).MarkerEdgeColor;
	
	
	% vertical line
	hpv(i) = plot(MergedData(i).Damping_1ocycles_fit_out.b*[1 1], ha.YLim);
	hpv(i).Color = hpf(i).Color;
	hpv(i).LineStyle = hpf(i).LineStyle;
	hpv(i).LineWidth = hpf(i).LineWidth;

    if flag_error
        dummy = confint(MergedData(i).Damping_1ocycles_fit_out);
        
        hplv(i) = plot(dummy(1,2)*[1 1], [0 20*ha.YLim(2)]);
        hplv(i).Color = hpf(i).Color;
        hplv(i).LineStyle = '--';
        hplv(i).LineWidth = hpf(i).LineWidth;
        
        hprv(i) = plot(dummy(2,2)*[1 1], [0 20*ha.YLim(2)]);
        hprv(i).Color = hpf(i).Color;
        hprv(i).LineStyle = '--';
        hprv(i).LineWidth = hpf(i).LineWidth;
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
ha.YLabel.String = '1/\tau_c, [cycles^{-1}]';
ha.XLabel.FontSize = 16;
ha.YLabel.FontSize = 16;


