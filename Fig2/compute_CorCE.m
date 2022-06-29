function [CorCE, theo_CorCE] = compute_CorCE(phi,binno,res_spks,mean_spks)

% phi, res_spks, mean_spks, binno


clear CorCE theo_CorCE

mean_mat_temp = nanmean(mean_spks,1);

CovCE = nancov(res_spks,'pairwise');
vartemp = nanvar(res_spks);

for i = 1:binno
    CovCE(i,i) = vartemp(i) - phi*mean_mat_temp(i);
    if (CovCE(i,i) < 0)
        CovCE(i,i) = nan;
    end
end

CorCE = nan(binno);
for i = 1:binno
    for j = 1:binno
        theo_CorCE(i,j) = sqrt(min(i,j)/max(i,j));
        CorCE(i,j) = CovCE(i,j)/sqrt(CovCE(i,i)*CovCE(j,j));
        if CorCE(i,j) > 1
            CorCE(i,j) = nan;
        end
    end
end


        
        
        
