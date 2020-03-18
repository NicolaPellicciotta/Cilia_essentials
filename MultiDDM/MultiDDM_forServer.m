addpath(genpath('D:\np451\'))
store_dir= 'S:\np451\odiri3\';
res_dir='D:\np451\odiri3_DDM\';

cd(store_dir);
d=dir('*CF1*.movie');
cd(res_dir)
for jj=1:numel(d);
    filename= strcat(store_dir,d(jj).name);
    %    if d(jj).bytes*1e-9 > 3% < 3 
    copyfile(filename,res_dir); 
    savename = strcat(d(jj).name(1:end-6),'.mat');
    disp(jj)
    boxsizes_vector = [16 32 64 128 256 512 1024];
    cilia = DDM_Analysis(filename);
    cilia.set_temperature(1);
    cilia.N_couple_frames_to_average = 200;
    cilia.VariableBoxSize_Analysis(boxsizes_vector);
    cilia.SAVAlike_CBF_measurement;
    cilia.gather_results;
    save(savename, 'cilia');
    clearvars -global
    clear cilia savename
    
    delete(d(jj).name)
 %       end
end