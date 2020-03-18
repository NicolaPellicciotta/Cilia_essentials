function [ Box ] = run_DDM_on_boxes( filename, BoxSize )
%run_DDM_on_boxes divides a Temika movie into boxes of size BoxSize (in px), runs
%DDM on each one and fits the Iqtau matrices
%   Detailed explanation goes here

% filename = '40X_beating.10Mar2015_09.17.47.movie'; %debugging
% BoxSize = 1028;                                    %debugging

%% input check

% if isempty(dir(filename)), error('File not found. You might be looking in the wrong folder'); end

if mod(BoxSize,2) %if odd BoxSize reduce of 1
    warning('BoxSize can''t be odd, reducing it of a unit.');
    BoxSize = BoxSize-1;
end

if BoxSize == 0
    error('BoxSize can''t be 0.');
end

%% open and load the movie file

fprintf('Opening the .movie file... ');

mo = moviereader(filename);

mo.width            = double(mo.width);             %uint32 to double cast to avoid funny integer maths
mo.height           = double(mo.height);            %uint32 to double cast to avoid funny integer maths

LastFrameToRead = mo.NumberOfFrames;

movie_correctly_read = false;   %control variable for failsafe movie reading
flag_movie_corrupted = false;   %control variable for failsafe movie reading

while movie_correctly_read == false %try to read the movie until it works
    try
        fs = mo.read([1 LastFrameToRead]);  %try to read the movie
        movie_correctly_read = true;        %if it works put flag to 1
    catch                                   %if it doesn't work, then
        LastFrameToRead = LastFrameToRead - 1;  %try to read again, ditching the last frame (I guess I could use a "finding" algorithm here)
        flag_movie_corrupted = true;        %set the flag that the movie was corrupted
    end                                     %and try again
    
end

fprintf('Done!');

if flag_movie_corrupted,    %warn user that the movie was corrupted
    fprintf('\n');
    warning(['Movie file was corrupted, only ',num2str(LastFrameToRead),' frames out of ',num2str(mo.NumberOfFrames),' were read. ']);
end

%% control on BoxSize and size of video: BoxSize can at most be as big as mo.width or mo.height

if BoxSize > mo.width || BoxSize > mo.height
    BoxSize = min(mo.width, mo.height);
    warning(['Impossible to run the analysis with the desired BoxSize, as it is bigger than at least one of the dimensions of the field of view. The analysis will be run on the maximum possible BoxSize, i.e. BoxSize = ',num2str(BoxSize)]);
end

%% prediction on sizes

N_Boxes_row = floor(mo.height / BoxSize);   %number of boxes that fit (vertically) in the field of view given BoxSize
N_Boxes_col = floor(mo.width / BoxSize);    %number of boxes that fit (horizontally) in the field of view given BoxSize
N_Boxes = N_Boxes_row * N_Boxes_col;        %number of boxes that fit (in total) in the field of view given BoxSize
N_Frames = size(fs,3);

max_mode = BoxSize/2;           %number of rows of Iqtau
max_tau = floor(N_Frames/2);    %number of cols of Iqtau

max_mode_fitted = floor(max_mode/2);    %maximum mode (higher q) fitted
max_tau_fitted = floor(max_tau);      %maximum lag time fitted.


%% initialisation of results structure

fprintf('\nInitialising data structure... ');

Box = struct(...
    'std_fs',    zeros(BoxSize,'single'),...
    'Iqtau',     zeros(max_mode, max_tau, 'single'),...
    'Frequency', zeros(max_mode_fitted,1),...
    'Amplitude', zeros(max_mode_fitted,1),...
    'Damping',   zeros(max_mode_fitted,1),...
    'Offset',    zeros(max_mode_fitted,1),...
    'GOF',       zeros(max_mode_fitted,1));
Box = repmat(Box, [N_Boxes, 1]);

fprintf('Done!');


%% choosing offsets to maximise output
% all the boxes touch each other, but if the sizes of the video are not
% multiples of BoxSize there is a bit of wiggle room to play with to
% maximise the signal we're analysing. this section just uses the standard
% deviation of the video to detect the most active region

fprintf('\nChoosing ROI placement... ');

std_fs = zeros(mo.height,mo.width,'single');

hpool = gcp;

parfor i=1:mo.height    %doing the std of the entire video is too RAM expensive, so for loop on the rows
    std_fs(i,:) = squeeze( std( single( fs(i,:,:) ) ,1,3) );
end

col_sum_std_fs = sum(std_fs,1); %sum of each column of std (vertically), it's a row vector
row_sum_std_fs = sum(std_fs,2); %sum of each row of std (horizontally), it's a column vector

col_span = N_Boxes_col * BoxSize;     %width of total DDM-analysed region (multiple of BoxSize)
row_span = N_Boxes_row * BoxSize;    %height of total DDM-analysed region (multiple of BoxSize)

%find where to place (in the horizontal direction) the [row_span, col_span] region to be DDM-analysed
for i = mo.width - col_span + 1 : -1 : 1  %sneaky allocation
    dummy(i) = sum(col_sum_std_fs(i:i+col_span-1));
end
[~,col_offset] = max(dummy);
clear dummy

%find where to place (in the vertical direction) the [row_span, col_span] region to be DDM-analysed
for i = mo.height - row_span + 1 : -1 : 1  %sneaky allocation
    dummy(i) = sum(row_sum_std_fs(i:i+row_span-1));
end
[~,row_offset] = max(dummy);
clear dummy

fprintf('Done!');


%% saving in each Box the relative bit of std_fs

fprintf('\nSaving portion of standard deviation of pixel intensity in time relative to each Box... ');
for i = 1:N_Boxes_row
    for j = 1:N_Boxes_col
        
        ii = sub2ind([N_Boxes_row, N_Boxes_col], i, j);
        Box(ii).std_fs = std_fs(row_offset + (i-1)*BoxSize : row_offset + (i)*BoxSize -1 ,...
                                col_offset + (j-1)*BoxSize : col_offset + (j)*BoxSize -1);
    end
end
fprintf(' Done!');


%% running DDM on each region of the video, storing Iqtaus in relative Box

fprintf('\nRunning DDM algorithm on each Box... ');
fprintf('\nAnalysing Box ');
for j = 1:N_Boxes_col
    for i = 1:N_Boxes_row
        
        ii = sub2ind([N_Boxes_row, N_Boxes_col], i, j);
        fprintf('%.6d / %.6d',ii,N_Boxes);
        Box(ii).Iqtau = DDM_core(fs(row_offset + (i-1)*BoxSize : row_offset + (i)*BoxSize -1 ,...
                                    col_offset + (j-1)*BoxSize : col_offset + (j)*BoxSize -1,:));
        fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b')
    end
end

fprintf(' Done!\n');


%% running fit_Iqtau_ourcilia on each Box.Iqtau

fprintf('\nFitting Iqtau matrix for each Box... ');
fprintf('\nFitting Box ');
for ii = 1:N_Boxes
        
    fprintf('%.6d / %.6d',ii,N_Boxes);
    [Box(ii).Frequency, Box(ii).Damping, Box(ii).Amplitude, Box(ii).GOF] = fit_Iqtau_ourcilia(Box(ii).Iqtau, max_mode_fitted, max_tau_fitted);
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b')
end

fprintf(' Done!\n');


