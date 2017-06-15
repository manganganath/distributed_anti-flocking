function gap = get_gap(r)

gap = zeros(length(r),length(r),2);

for a=2:1:length(r)
    for b=1:1:a-1
        gap(a,b,:)=r(b,:)-r(a,:);
    end 
end

gap(:,:,1) = gap(:,:,1) - gap(:,:,1)';
gap(:,:,2) = gap(:,:,2) - gap(:,:,2)';

return