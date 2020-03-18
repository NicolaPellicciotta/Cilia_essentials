function res = AutoCorr(cc,ll)

if size(cc,1) < (ll-1)
    disp('ciccia...')
    return
else
    if length(cc) < 2*ll
        disp('attenzione! pochi dati ...');
    end
end

MCC2 = mean(cc).^2; %row vector

for i = ll:-1:1
   res(i,:) = mean(cc(1:(end+1-i),:).*cc(i:end,:)) - MCC2;  % each of these is a row
end

% res = res./res(1,:); % thanks to implicit expansion, each column is normalised
res = res ./ var(cc);