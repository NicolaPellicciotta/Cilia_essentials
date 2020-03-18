function [bs,D,eD,dD]=Sigmoid(cilia)

Nbox= numel(cilia.Results);
px2mu=0.14;
kk=1;
for b=1:Nbox
    bsz= cilia.Results(b).BoxSize;
    q_limits_1oum = [1.95 2.75];
%    q_limits_1oum = [4 10];
    mod_lim= floor(q_limits_1oum*px2mu* bsz /(2*pi)); 
    if all(mod_lim)==1;
    
    ind_good= cilia.Results(b).ind_good_boxes(:);
    Damping=zeros([1,sum(ind_good)]);
    cc=1;
    for j=1:numel(cilia.Results(b).Box);
        if ind_good(j)==1;
        Damping(cc) = nanmedian(cilia.Results(b).Box(j).Damping(mod_lim(1):mod_lim(2))); 
        cc=cc+1;
        end
    end
    
    q_limitis_1opx= q_limits_1oum*px2mu;
    
    D(kk)= nanmedian(Damping);
    eD(kk)= nanmedian( abs(Damping -D(kk)))./sqrt(numel(Damping));
    dD{kk}= Damping; 
    bs(kk)=cilia.Results(b).BoxSize;
%    histogram(dD{1});
    kk=kk+1;
    end   
   
end
%figure(2)
% plot(bs,D,'-o');