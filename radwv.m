function meanrad = radwv(refl_full,spec,wv)

if size(wv,1) == 1     % narrow band
    meanrad = nan(size(wv));

    for i = 1 : length(wv)
        meanrad(i) = mean(refl_full(spec>=wv(i) & spec<wv(i)+0.1));
    end
    
elseif size(wv,1)>1 && size(wv,2)==2   % wide band
    meanrad = nan([1,size(wv,1)]);
    
    for i = 1 : size(wv,1)
        meanrad(i) = mean(refl_full(spec>=wv(i,1) & spec<=wv(i,2)));
    end
end
    