function [F,BW] = remove_debris(F_temp,area_min,f_lim)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
f=ceil(F_temp);
[G,idfreq]=findgroups(f(:));
BW=zeros(size(f));
for jj=1:numel(idfreq)
    if idfreq(jj)>=f_lim(1) & idfreq(jj)<f_lim(2)
    bw=zeros(size(f));
    bw(f==idfreq(jj))=1;bw=logical(bw);
    bw2=bwareaopen(bw,area_min);
    BW(bw2)=1;
    end
end
BW=logical(BW);
F=F_temp;
F(~BW)=nan;
end

