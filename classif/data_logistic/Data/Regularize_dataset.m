% Regularize dataset

datasetname = 'sido0_data';


load(datasetname)

for feature_idx = 1:nFeatures
    display([num2str(round(100*feature_idx/nFeatures)), '%'])
    vec = X(feature_idx,:);
    nnzvec = find(vec);
    maxvec = max(abs(vec));
    if(maxvec ~= 0)
        X(feature_idx,nnzvec) = X(feature_idx,nnzvec)/maxvec;
    end
    
end

clear feature_idx max_feature vec maxvec nnzvec

save([datasetname, '_normalized']);