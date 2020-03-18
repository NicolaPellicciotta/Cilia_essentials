clear all;
close all;

path = pwd;
list=dir('alg*.movie');

for i=1:numel(list)
    
filename = list(i).name

% dir_string = '*.tiff';
% frame_stack = load_tiff_sequence(dir_string, path);

tic;
movie_object = moviereader(filename);
video_struct = movie_object.read;
timestamp=[video_struct(:).relative_timestamp_sec];
video_struct = [];
frames2s = mean(diff(timestamp));

toc;
Iqtau = DDM_core_GPU(movie_object);
toc;
parsave([ filename(1:end-5),'mat'],Iqtau,timestamp,frames2s);
toc;


end
