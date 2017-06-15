function cell_map = fuse_record(cell_map,nbr)

for i=2:1:size(nbr,1)
    for j=1:1:i-1
        if nbr(i,j)==1
            max_time = max(cell_map(:,:,2*i-1),cell_map(:,:,2*j-1));
            temp_cell_i = cell_map(:,:,2*i);
            temp_cell_i(cell_map(:,:,2*i-1) ~= max_time) = j;
            cell_map(:,:,2*i-1) = max_time;
            cell_map(:,:,2*j-1) = max_time;
            cell_map(:,:,2*i) = temp_cell_i;
            cell_map(:,:,2*j) = temp_cell_i;
        end
    end
end

return