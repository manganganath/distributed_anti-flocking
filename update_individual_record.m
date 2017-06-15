function cell_map = update_individual_record(cell_map,map_pos,x,T,r_s)

for r=1:1:size(x,1)
    other_dist = sqrt((x(r,1)-map_pos(:,:,1)).^2+(x(r,2)-map_pos(:,:,2)).^2);
    temp_map = cell_map(:,:,2*r-1);
    temp_ind = cell_map(:,:,2*r);
    temp_map(other_dist<=r_s)=T;
    temp_ind(other_dist<=r_s)=r;
    cell_map(:,:,2*r-1)=temp_map;
    cell_map(:,:,2*r)=temp_ind;
end