% define angle map 
jj = repmat([1:frame_size],frame_size,1);
ii=jj';
cc = max_q+1;



angle_map = atand( (ii-cc)./(jj-cc) );

%pick only the angles that we care about
for i=1:1024
    for j=1:1024
        value = angle_map(i,j);
        if (value > -10) & (value < 10)
            angle_map(i,j) = 2;
        elseif (value > 80) | (value < -80)
            angle_map(i,j) = 1;
        else
            angle_map(i,j) = 3;
        end
    end
end



imagesc(angle_map); colorbar; shg

angle_map = fftshift(angle_map);