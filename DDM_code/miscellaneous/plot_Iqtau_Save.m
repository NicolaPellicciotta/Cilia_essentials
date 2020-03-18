function [ ] = plot_Iqtau_Save(px2mum, min_q_fit, max_q_fit, index, foo, foo2, fitdir)
%% Saves the figure labelled by 'index' under 'fitdir'

%%
% Which fit?
switch index
    case 1
        nme = 'tauq_';
        Conv = px2mum^2;
    case 3
        nme = 'omega_';
        Conv = px2mum;
end

% Save under
% userdir/DDM_core2_[Th]_[mmdd_HHMM]/fit_[max_q]_[max_tau]_[func]_[mmdd_HHMM]/[nme]_[min_q_fit]-[max_q_fit]_[alpha]_[A]_[mmdd_HHMM].[ext]
t = datestr(now, 'mmdd_HHMM');
plotname = [nme,num2str(min_q_fit),'-',num2str(max_q_fit), ...
                '_',num2str(foo2.alpha),'_',num2str(foo2.A*Conv),'_',t];
plotext = 'fig';
dataext = 'mat';
plotfilename = strjoin({plotname,plotext},'.');
datafilename = strjoin({plotname,dataext},'.');
plotpath = strjoin({fitdir,plotfilename}, filesep);
datapath = strjoin({fitdir,datafilename}, filesep);
h = figure(index);
variables = {'foo','foo2','min_q_fit','max_q_fit','px2mum'};  %,'tau_q_array','q','tau'};

savefig(h, plotpath, 'compact');
disp(['Saved figure (',num2str(index),').'])

save(datapath,variables{:});
disp(['Saved fit data for figure (',num2str(index),').'])

end