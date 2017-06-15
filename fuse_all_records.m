function fused_scan_record = fuse_all_records(cell_map,num_agents,fused_scan_record)

max_time = max(fused_scan_record(:,:,1),cell_map(:,:,1));
temp_fused = fused_scan_record(:,:,2);
temp_fused(max_time~=fused_scan_record(:,:,1)) = 1;
fused_scan_record(:,:,2) = temp_fused;
fused_scan_record(:,:,1) = max_time;

for i=3:2:num_agents*2
    max_time = max(fused_scan_record(:,:,1),cell_map(:,:,i));
    temp_fused = fused_scan_record(:,:,2);
    temp_fused(max_time~=fused_scan_record(:,:,1)) = (i+1)/2;
    fused_scan_record(:,:,2) = temp_fused;
    fused_scan_record(:,:,1) = max_time;
end