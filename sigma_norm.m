function z_n = sigma_norm(z,efs)

z_n = (1/efs)*(sqrt(1+efs*squd_norm(z))-1);
