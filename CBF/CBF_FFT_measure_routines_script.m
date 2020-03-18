    % set the range of frequency that you expect, this is [4,15] for airways
    f_lim=[4,15];  
    % only calculate frequency of cilia. Cilia are identified through a standard
    % deviation of pixels intensity.
    % First set the threeshold  probaility to be a cilia 
    %the smaller is threshold and the higher is the number of cilia that it fin

    threshold_cilia=0.5; %%% between 0 and 1, play a bit with this
    box_size=4; % fft will be averaged on a box_size
    area_min=10; % mean area of a ciliated cell in pixels
    [F_temp,s,BW_temp] = find_cilia(mo,f_lim,threshold_cilia,box_size);
    [F4,BW] = remove_debris(F_temp,area_min,f_lim);
    ss=imresize(BW,box_size,'nearest');    
    IF= imresize(F4,box_size,'nearest');