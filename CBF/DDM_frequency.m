%% use DDM to find beating frequency. this is very slow, not reccommended. 

filename= '40X_RR.12Dec2018_12.20.53_M_T25_A3.movie'
box_size=32;
cilia = DDM_Analysis_nico(filename);
cilia.VariableBoxSize_Analysis(box_size);
save([cilia.Filename(1:end-5),'mat'],'cilia');  %%% save result

cilia.gather_results;

boxes=[];
for bb= 1:numel(cilia.Results); boxes(bb)= cilia.Results(bb).BoxSize;
end
bsi = find(boxes == box_size);
freq =cilia.Results(bsi).MedianFrequencyVec