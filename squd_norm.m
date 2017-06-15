function s_n = squd_norm(z)

if size(z,3)==1
    s_n = z^2;
else
    s_n = z(:,:,1).^2 + z(:,:,2).^2;
end