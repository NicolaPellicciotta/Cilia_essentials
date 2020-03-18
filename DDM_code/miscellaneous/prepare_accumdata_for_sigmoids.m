%% D:\Data\Cilia\Analysis\2017_08_04--07

cch

cd('D:\Data\Cilia\Analysis\2017_08_04--07');



% prepare variables for populating data structure
%
analysis_folder = {'D:\Data\Cilia\Analysis\2017_08_04--07\2017_08_04_CF_LUMIIVA_00h_Analysis';...
    'D:\Data\Cilia\Analysis\2017_08_04--07\2017_08_05_CF_LUMIIVA_24h_Analysis';...
    'D:\Data\Cilia\Analysis\2017_08_04--07\2017_08_06_CF_LUMIIVA_48h_Analysis';...
    'D:\Data\Cilia\Analysis\2017_08_04--07\2017_08_07_CF_LUMIIVA_72h_Analysis'};
figures_folder = 'Figures_sigmoids';
imaging_string = '40X_BF_CF';

sampletypes = {'LumacaFTOR';
    'Orkambi';...
    'T1a';...
    'Control'};

timepoints = {'00h';...
    '24h';...
    '48h';...
    '72h'};

donors = {'d1';'d2';'d3'};

inserts = {''};

positions = '';

% use default quick one, which is in this case the correct one
% boxsizes_vector = [16 32 48 64 96 128 160 192 224 256 340 512 1024];
boxsizes_vector = [32 64 128 256 512 1024];

% use the default one, works with this magnification/camera combo
q_limits_1oum = [];

flag_dryrun = false;
flag_recalculate_goodboxes = true;

[SampleType] = populate_DDM_Sigmoids_struct( analysis_folder, figures_folder,...
    imaging_string, sampletypes, timepoints, donors, inserts, positions,  boxsizes_vector, q_limits_1oum,...
    flag_dryrun, flag_recalculate_goodboxes);

% return
% save
save 'AccumData_from_func_newboxes_nosparse.mat';
mkdir(figures_folder)