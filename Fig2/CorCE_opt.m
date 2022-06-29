function err = CorCE_opt(phi,binno,res_spk,mean_spk)

[CorCE, theo_CorCE] = compute_CorCE(phi,binno,res_spk,mean_spk);

%% fisher z-transformation of r to z

for i = 1:length(CorCE)
    for j = 1:length(CorCE)        
        if i==j
            z_CorCE(i,j) = 1;
            z_theo_CorCE(i,j)  = 1;
        else
            z_CorCE(i,j) = 0.5*(log(1+CorCE(i,j))-log(1-CorCE(i,j)));
            z_theo_CorCE(i,j) = 0.5*(log(1+theo_CorCE(i,j))-log(1-theo_CorCE(i,j)));
        end
    end
end

err = sum(sum((z_CorCE-z_theo_CorCE).^2));
fprintf('phi=%.2f\t err=%g\n',phi,err);