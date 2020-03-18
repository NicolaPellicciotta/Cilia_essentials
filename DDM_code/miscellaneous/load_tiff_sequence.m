function [frame_stack, frame_size, frame_class] = load_tiff_sequence ( dir_string, path, N_frames)

%% input controls
if nargin <3
    N_frames = [];
    if nargin<2
        path = pwd;
    end
end

%% detects the names of the tiff series

frame_names = dir(fullfile(path,dir_string));
if isempty(frame_names)
    disp('Error: no tiffs in current directory');
    return;
end

if isempty(N_frames),
    N_frames = size(frame_names,1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_frames = 1024;     %for quicker debugging, I just use a few frames
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% loads and extracts metadata from the first image

frame = imread(fullfile(path,frame_names(1).name));
frame_info = whos('frame');
frame_size = frame_info.size;
frame_size = min(frame_size);
frame_class = frame_info.class;
clear frame

%% loads the tiffs into a width x height x number_of_frames matrix

frame_stack = zeros([frame_size, frame_size, N_frames], frame_class);
for i=1:N_frames
    frame_stack(:,:,i) = imread(fullfile(path,frame_names(i).name), 'PixelRegion', {[1 frame_size], [1 frame_size]});
end

frame_stack = gpuArray(frame_stack);
end


