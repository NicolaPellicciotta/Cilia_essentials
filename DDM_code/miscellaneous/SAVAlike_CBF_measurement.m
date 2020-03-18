function [frequency_map, median_frequency, std_video,pks_map] = SAVAlike_CBF_measurement(filename, flag_ROI)

if nargin < 1
    filename = uigetfile('*.movie');
end

if nargin < 2
    flag_ROI = false;
end
% filename = 'D:\Data\Cilia\Data\2015_02_20\cilia39.20Feb2015_12.10.32.movie';

mo = moviereader(filename);
fs = mo.read;

bsz = 4; %size of the binning box

%% select ROI
hpool = gcp;

std_fs = zeros(mo.height, mo.width, 'single');

parfor i=1:mo.height
    std_fs(i,:) = squeeze(std(single(fs(i,:,:)),1,3));
end

if flag_ROI
    figure;
    imagesc(std_fs); colormap jet;
    
    rect = round(getrect(gcf));
        
    rect(3) = bsz*ceil(rect(3)/bsz);
    rect(4) = bsz*ceil(rect(4)/bsz);

    fs = fs(rect(2):rect(2)+rect(4)-1, rect(1):rect(1)+rect(3)-1,:);
end

[height, width, N_frames] = size(fs);

FrameRate = mo.FrameRate;
%%
binningmap = col2im(repmat(1:floor(width*height/bsz^2),bsz*bsz,1), bsz.*[1 1], bsz*[floor(height/bsz) floor(width/bsz)],'distinct');

binnedfs = zeros(floor( [height/bsz width/bsz N_frames]), 'single');

binna = @(i) reshape( 1./bsz^2.*accumarray(binningmap(:),single(reshape(fs(:,:,i),width*height,1))), [floor(height/bsz), floor(width/bsz)]);
for i=1:N_frames
    binnedfs(:,:,i) = binna(i);
end



% clear fs

bheight = floor(size(binningmap,1)/bsz);
bwidth  = floor(size(binningmap,2)/bsz);

frequency_map = zeros(bheight, bwidth,'double');

figure
tic;

binnedfs = reshape(binnedfs,[bheight*bwidth N_frames]);
fy = abs(fft(binnedfs,[],2));
spectrum1 = fy(:,1:ceil(N_frames/2)); %next three line from algorithm to detect peaks in signal with harmonics
spectrum2 = fy(:,1:2:N_frames);
product_spectrum = spectrum1.*spectrum2;

pks_map = zeros(size(frequency_map),'double');
parfor ii=1:bheight*bwidth
    
    [pks,loc,~,proms] = findpeaks( double(product_spectrum(ii,:)) );
%     [~,loc,~,proms] = findpeaks( double(spectrum1(ii,:)) ); %good only for very clean spectra
    [pkmax,ind] = max(pks);
    %     [loc, pks] = peakseek(product_spectrum(ii,:));
    %     [~,ind] = max(pks);
    frequency_map(ii) = (loc(ind)+1)/(N_frames).* FrameRate; %converting to Hz (? to be checked)
    pks_map(ii) = pkmax;
end

hist(frequency_map(:),100);
median_frequency = median(frequency_map(:));
std_video = std_fs;
toc
end