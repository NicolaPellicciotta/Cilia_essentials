%% Multi-DDM analysis by L.Feriani but as interpreted by N.Pellicciotta.  
% the first way is to just use the GUI that Luigi made for biologist
Analyse_Epithelix_func
% here you can run Multi-DDM on as amny video as you want.
% the Alternative version is doing this:

store_dir= 'S:\np451\odiri3\';  %%% folder where the data are
res_dir='D:\np451\odiri3_DDM\'; %%% folder where you want to save your results

cd(store_dir);
d=dir('*.movie');
for jj=1:numel(d);
    filename= strcat(store_dir,d(jj).name);
    disp(jj)
    boxsizes_vector = [16 32 64 128 256 512 1024];
    cilia = DDM_Analysis(filename);
    cilia.set_temperature(1);
    cilia.N_couple_frames_to_average = 200;
    cilia.VariableBoxSize_Analysis(boxsizes_vector);
    cilia.SAVAlike_CBF_measurement;
    cilia.gather_results;
    cd(res_dir)
    save(savename, 'cilia');
    clearvars -global
    clear cilia savename
    cd(store_dir)
 
end

%% After you have your results, you need to join all the info from the same samples 
% in the variable SampleType, then saved in AccummData.mat
% this works only if you wrote the name of your results .mat in such a way
% that the string corresponding to drug treatment or temperature is
% coherent. 

analysis_folder = 'D:\np451\odiri2_DDM'
figures_folder = 'Figures_sigmoids';
imaging_string = '40X';

%%%% the best way is to vary thighs at one, so decide a sampletype and vary the timespoints (the two are indipendednt) 


sampletypes = {'_D','_L','_R'}   %%%% divide for drug treatment

timepoints = {'0h','48h'};  %% or can be temperature
donors={'_CF1_';'_CF2_';'_CF3_'};  %%% different donors
%inserts = {'1_','2_','3_'};
%positions = {''};
%donors={''};
inserts = {''};
positions = {''};


% use default quick one, which is in this case the correct one
 boxsizes_vector = [16 32 64 128 256 512 1024];
%boxsizes_vector = [16 32 64 128 256 512 1024];

% use the default one, works with this magnification/camera combo
q_limits_1oum = [];

flag_dryrun = false;
flag_recalculate_goodboxes = true;

[SampleType] = populate_DDM_Sigmoids_struct( analysis_folder, figures_folder,...
    imaging_string, sampletypes, timepoints, donors, inserts, positions,  boxsizes_vector, q_limits_1oum,...
    flag_dryrun, flag_recalculate_goodboxes);

save 'AccumData.mat';
mkdir(figures_folder)

% now for plots see the MultiDDM_sigmoids.m script
