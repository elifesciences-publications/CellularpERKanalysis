function [ idx_subset, new_idx ] = resampleindices( fulllengthdata, fulllengthindices, nplot, pc)

    if nplot == 'auto',
        
    nplot = (pc/100) * length(fulllengthdata);
    
    else
    end

     perm = randperm(length(fulllengthdata)); 
     idx_subset = perm(1:nplot); %get indexes
     Lia = ismember(idx_subset, fulllengthindices);
     new_idx = idx_subset(Lia.*idx_subset > 0);

end

