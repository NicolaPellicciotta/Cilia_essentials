function [ binning_map ] = create_binning_map( output_size, box_size )
%create_binning_map Makes stuff with indices to get you a tiled matrix with
%different numbers for each tile
%   this will work with 2016b, not sure earlier versions

%% input check

if isscalar(output_size)
    output_size = repmat(output_size,1,2);
end

if isscalar(box_size)
    box_size = repmat(box_size,1,2);
end

if any(mod(output_size,box_size))
    error('output_size has to be a multiple of box_size')
end

%% using implicit indices expansion 

binning_map = ceil((1:output_size(1))/box_size(1))' +...
    output_size(1)/box_size(1) * (ceil((1:output_size(2))/box_size(2))-1);


end

