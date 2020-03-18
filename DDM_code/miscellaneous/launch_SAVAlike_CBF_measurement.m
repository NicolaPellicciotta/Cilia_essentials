clear
close all

datadirname     = 'D:\Data\Cilia\Data\2015_02_20';
analysisdirname = 'D:\Data\Cilia\Analysis\2015_02_20';

filenames = dir(fullfile(datadirname,'*.movie'));

for i=8:numel(filenames)
    display(filenames(i).name)
    try
        frequency_map = SAVAlike_CBF_measurement(fullfile(datadirname,filenames(i).name));
    catch
        continue
    end
    try
        save(fullfile(analysisdirname,['SAVAlike_CBF_',filenames(i).name(1:end-5)]),'frequency_map');
    catch
        mkdir(analysisdirname);
        save(fullfile(analysisdirname,['SAVAlike_CBF_',filenames(i).name(1:end-5)]),'frequency_map');
    end
end
