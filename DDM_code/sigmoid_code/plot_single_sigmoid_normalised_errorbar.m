function [] = plot_single_sigmoid_normalised_errorbar(MergedData, flag_confint)
%plot_single_sigmoid_normalised_errorbar takes a MergedData structure (as output from
%merge_SampleType_data) and plots a single sigmoid with the fitting line
%and errorbar

if nargin < 2 || isempty(flag_confint)
    flag_confint = false;
end

% prepare figure
hf = figure;

% prepare axes
ha = axes;
box on;
hold on;

% axes limits
setsemilogx
ha.YLim = [0 3];
ha.XLim = [1e0 1e5];
ha.XTick = logspace(0,5,6);

% errorbar
he = errorbar(MergedData.window_area_um2, MergedData.med_Damping_1ocycles,...
    MergedData.ler_Damping_1ocycles, MergedData.uer_Damping_1ocycles);
he.LineWidth = 1.2;
he.LineStyle = 'none';
he.Marker = 'o';
he.MarkerFaceColor = 'w';

% fit
xx = logspace(0,5,1e3);
hpf = plot(xx, MergedData.Damping_1ocycles_fit_out(xx));
hpf.LineWidth = 1.2;

% vertical line - is now the shoulder
hpv = plot(exp(2).*MergedData.Damping_1ocycles_fit_out.b*[1 1], ha.YLim);
hpv.Color = hpf.Color;
hpv.LineStyle = hpf.LineStyle;
hpv.LineWidth = hpf.LineWidth;

if flag_confint
    dummy = confint(MergedData.Damping_1ocycles_fit_out);
    dummy = dummy(:,2);
    dummy = dummy .* exp(2);
    xp = dummy([1 2 2 1]);
    yp = [0 0 1e6 1e6];
    hpp = patch(xp, yp, 'red',...
        'FaceColor', hpv.Color,...
        'FaceAlpha', 0.2,...
        'EdgeColor', 'none');
end

 
% labels and titles
ha.XLabel.String = 'DDM Window Area, [\mum^2]';
ha.YLabel.String = '1/\tau_c, [cycles^{-1}]';
ha.XLabel.FontSize = 16;
ha.YLabel.FontSize = 16;


