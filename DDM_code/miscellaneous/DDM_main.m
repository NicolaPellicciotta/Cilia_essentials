clear all;
close all;

path = pwd;
list=dir('sample*.movie');

for i=1:numel(list)
    
filename = list(i).name

% dir_string = '*.tiff';
% frame_stack = load_tiff_sequence(dir_string, path);

tic;
movie_object = moviereader(filename);
video_struct = movie_object.read;
frame_stack = cat(3,video_struct(:).IM);
timestamp=[video_struct(:).relative_timestamp_sec];
frames2s = mean(diff(timestamp));
video_struct = [];

toc;
Iqtau = DDM_core(frame_stack,false);
toc;
save([ filename(1:end-5),'mat'],'Iqtau','timestamp','frames2s');
%fit_Iqtau_colloids;
toc;

clearvars -except path list

end