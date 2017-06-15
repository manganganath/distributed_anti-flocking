%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% An implementation of
% @article{ganganath2016anti-flocking,
% title={Distributed Anti-Flocking Algorithms for Dynamic Coverage of Mobile Sensor Networks},
% author={Nuwan Ganganath and Chi-Tsun Cheng and Chi K. Tse},
% journal={IEEE Transactions on Industrial Informatics},
% year={2016},
% volume={12},
% number={5},
% pages={1795--1805},
% publisher={IEEE}}
%
% Tested on MATLAB R2011b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clear all; clc;

% Parameters
num_agents = 5; % Number of agents (min = 3)
t_gap = 0.1; % time gap between two iterations (s)
r_c = 20;  % communication radius
r_s = 5; % sensing radius
v_0 = 1; %initial velocity
efs = 0.1;
c1 = 10;  %decentering
c2 = 0.5;  %actraction to the target
c3 = 0.8;  %velocity feedback
rho = 0.2;
sigma1 = 0.04;
sigma2 = 0.01;
d_a = 1.9*r_s;
d_o = 0.9*r_s;
map_width = 50; %width of a squre map
map_res = 0.5; %width of a grid

% Defining vectors
x = zeros(num_agents,2);        %current position
x_1 = zeros(num_agents,2);      %previous position
v_1 = zeros(num_agents,2);      %previous velocity
tra_dist = zeros(num_agents,1); %travel distance

cell_map = zeros(floor(map_width/map_res),floor(map_width/map_res),num_agents*2); %cell_map structure

%initializing the grid centers
map_pos = zeros(size(cell_map,1),size(cell_map,2),2);
for g_i=1:1:size(cell_map,1)
    for g_j=1:1:size(cell_map,2)
        map_pos(g_i,g_j,1)= map_res*(g_i-1/2);
        map_pos(g_i,g_j,2)= map_res*(g_j-1/2);
    end
end
 
%initializing obstacles
obs_map = zeros(size(cell_map,1),size(cell_map,2));
obs_map(65:80,15:30)=-Inf;
obs_map(sqrt((map_pos(:,:,1)-10).^2+(map_pos(:,:,2)-10).^2)<3)=-Inf;
obs_map(sqrt((map_pos(:,:,1)-30).^2+(map_pos(:,:,2)-40).^2)<5)=-Inf;

%initializing agents
map_pos_x = map_pos(:,:,1);
map_pos_y = map_pos(:,:,2);
map_free_x = map_pos_x(obs_map~=-Inf); 
map_free_y = map_pos_y(obs_map~=-Inf);
free_ind = datasample(1:1:length(map_free_x),num_agents);
x_1(:,:) = [map_free_x(free_ind) map_free_y(free_ind)];

x_r = x_1;                      %target position
x_r_t = zeros(num_agents,3);    %target last visited time
[I,J] = ind2sub(size(map_pos(:,:,1)),free_ind);
x_r_t(:,1:2) = [I',J'];

ini_v_t = 2*pi*rand(num_agents,1);
ini_v_s = v_0*rand(num_agents,1);
v_1(:,1) = ini_v_s.*cos(ini_v_t);
v_1(:,2) = ini_v_s.*sin(ini_v_t);

%defining agents for obstacles
x_obs = [map_pos_x(obs_map==-Inf) map_pos_y(obs_map==-Inf)];
num_obs_agents = size(x_obs,1);
obs_nbr = zeros(num_obs_agents,num_agents);

fused_scan_record = zeros(size(cell_map,1),size(cell_map,2),2);
fused_scan_record(:,:,1) = fused_scan_record(:,:,1) - obs_map;

counter = 1;

r_c_sigma = sigma_norm(r_c,efs);
r_s_sigma = sigma_norm(r_s,efs);

% main loop
while(true)    
    u_d = zeros(num_agents,2);      %control input for decentralizing 
    u_o = zeros(num_agents,2);      %control input for obstacle avoidance 
    
    % calculating adjacency matrix 
    dist_gap = get_gap(x_1);
    dist_2 = squd_norm(dist_gap);
    dist = sqrt(dist_2);
    adj = action_func(dist,d_a,c1);
    nbr = zeros(num_agents);
    nbr(dist<=r_c) = 1;
    nbr = nbr - diag(diag(nbr)); %set diagonal to 0
    adj = nbr.*adj;
    
    T=counter*t_gap;
        
    cell_map = update_individual_record(cell_map,map_pos,x_1,T,r_s);
    
    if sum(sum(nbr))>0
        cell_map = fuse_record(cell_map,nbr);
    end
    
    % fusing cell_map for calculations
    fused_scan_record = fuse_all_records(cell_map,num_agents,fused_scan_record);     
        
    rem_map=sum(sum(fused_scan_record(:,:,1)==0));
    cur_covg=sum(sum(fused_scan_record(:,:,1)==T));
    [map_w,map_h,~]=size(fused_scan_record);
    cumu_covg = (map_w*map_h-rem_map-num_obs_agents)/((map_w*map_h-num_obs_agents)*0.01); 
    inst_covg = cur_covg/((map_w*map_h-num_obs_agents)*0.01);       

    %Selfishness
    u_e = c2*(x_r-x_1)-c3*v_1;  
        
    %decentering
    for a=1:1:num_agents
        for b=1:1:num_agents
            u_d(a,:) = u_d(a,:) + adj(a,b)*(x_1(b,:)-x_1(a,:))/sqrt(efs+dist_2(a,b)); 
        end
    end
        
    %obstacle avoidance
    if(num_obs_agents>0)
        for i=1:1:num_agents
            obs_nbr(:,i)=sqrt((x_obs(:,1)-x_1(i,1)).^2+(x_obs(:,2)-x_1(i,2)).^2);
        end
        obs_nbr(obs_nbr<=d_o)=1;
        obs_nbr(obs_nbr>d_o)=0;
        
        for a=1:1:num_agents
            obs_nbr_ind = find(obs_nbr(:,a)==1);
            if ~isempty(obs_nbr_ind)
                for b=1:1:length(obs_nbr_ind)
                    o_gap = zeros(1,1,2);
                    o_gap(:,:,1)=x_obs(obs_nbr_ind(b),1)-x_1(a,1);
                    o_gap(:,:,2)=x_obs(obs_nbr_ind(b),2)-x_1(a,2);
                    o_dist_2 = squd_norm(o_gap);
                    u_o(a,:) = u_o(a,:) + action_func(sqrt(o_dist_2),d_o,c1)*...
                        (x_obs(obs_nbr_ind(b),:)-x_1(a,:))/sqrt(efs+o_dist_2); %obstacle avoidance
                end
                
                if acos(dot(u_o(a,:),v_1(a,:))/(norm(v_1(a,:))*norm(u_o(a,:))))<pi/2
                    u_o(a,:)=[0 0];    
                end                             
            end
        end        
    end

    u = u_d + u_o + u_e;
    
    %calculating the new position and velocity
    x(:,1:2) = x_1(:,1:2) + v_1*t_gap;
    v = v_1 + u*t_gap;
    
    v_1 = v;
    x_pre = x_1;
    x_1 = x;
    tra_dist = tra_dist + sqrt((x(:,1)-x_pre(:,1)).^2+(x(:,2)-x_pre(:,2)).^2);
    
    %calculate goal distances
    goal_dist = zeros(num_agents,num_agents);
    for s=1:1:num_agents %agent
        for t=1:1:num_agents %agent
            goal_dist(s,t) = sqrt((x(s,1)-x_r(t,1)).^2+(x(s,2)-x_r(t,2)).^2);
        end
    end

    for s=1:1:num_agents
        obs_map_temp = obs_map;    
        temp_cell_map = cell_map(:,:,2*s-1);     
        recalculate = 0;
        if goal_dist(s,s)<r_s  % if target is covered by assigned agent
            recalculate = 1;
        elseif temp_cell_map(x_r_t(s,1),x_r_t(s,2)) ~= x_r_t(s,3) % if target is covered by someone else
            recalculate = 1;
        else % if target is close to another target
            for i=1:1:num_agents
                if nbr(s,i)==1
                    if (sqrt((x_r(i,1)-x_r(s,1)).^2+(x_r(i,2)-x_r(s,2)).^2) < r_s) && (goal_dist(i,i) < goal_dist(s,s))
                        obs_map_temp(sqrt((map_pos(:,:,1)-x_r(i,1)).^2+(map_pos(:,:,2)-x_r(i,2)).^2)<r_s) = -Inf;
                        recalculate = 1;
                    elseif (goal_dist(i,s) < goal_dist(s,s)) && (goal_dist(s,i) < goal_dist(i,i)) && (recalculate == 0)
                        temp_x_r = x_r(i,:);
                        x_r(i,:) = x_r(s,:);
                        x_r(s,:) = temp_x_r;
                    end
                end
            end
        end

        if recalculate == 1               
            %calculate distance to each grid 
            tar_dist = sqrt((x(s,1)-map_pos(:,:,1)).^2+(x(s,2)-map_pos(:,:,2)).^2);
            %calculate minimum distance to each grid from neighbors
            nbr_tar_dist = Inf(size(tar_dist));
            for i=1:1:num_agents
                if nbr(s,i)==1
                    nbr_tar_dist = min(nbr_tar_dist,sqrt((x(i,1)-map_pos(:,:,1)).^2+(x(i,2)-map_pos(:,:,2)).^2));
                end
            end

            obs_map_temp(tar_dist~=min(nbr_tar_dist,tar_dist))=-Inf;
            pre_tar_dist = sqrt((x_r(s,1)-map_pos(:,:,1)).^2+(x_r(s,2)-map_pos(:,:,2)).^2);
            map_d = (T-temp_cell_map).*(rho+(1-rho)*exp(-sigma1*tar_dist-sigma2*pre_tar_dist)); 
            map_d1 = map_d.*(obs_map_temp+1);
            [~,ind] = max(map_d1(:));
            [I,J]=ind2sub([size(map_d,1) size(map_d,2)],datasample(ind,1));
            x_r(s,:) = map_pos(I,J,:);
            x_r_t(s,:) = [I,J,temp_cell_map(I,J)];
        end
    end
       
    %plotting       
    h1=figure(1);
    subplot(1,2,1)
    axis square;box on;
    cc=hsv(num_agents+1);

    for i=1:1:num_agents       
        f1 = findall(gca,'Color',cc(i+1,:));
        f2 = findall(gca,'FaceColor',cc(i+1,:));
        delete(f1);delete(f2);
    end              
    delete(findall(gca, 'color', 'k'));
    delete(findall(gca, 'color', 'r'));
    refresh;    


    gplot(nbr,x_1,'r-');
    max_abs = sqrt(v_1(:,1).^2+v_1(:,2).^2);
    x_e = x_1+2*v_1./[max_abs max_abs];       

    hold on
    for i=1:1:num_agents       
        plot(x_r(i,1),x_r(i,2),'h','Color',cc(i+1,:),'MarkerSize',6);
        plot(x_1(i,1),x_1(i,2),'o','Color',cc(i+1,:),'MarkerSize',5);
        if T>0.2
        plot([x_pre(i,1),x(i,1)],[x_pre(i,2),x(i,2)],'Color',abs(cc(i+1,:)-0.1));
        end
    end          

    xlabel('x (m)');ylabel('y (m)');
    xlim([-map_width/10 map_width+map_width/10]);ylim([-map_width/10 map_width+map_width/10])  
    title(['locations of \alpha-agents at t = ',num2str(T-0.2),'s']);
    if(num_obs_agents>0)
        plot(x_obs(:,1),x_obs(:,2),'ks','MarkerSize',5); 
    end  

    subplot(1,2,2)
    surf(map_res:map_res:map_width,map_res:map_res:map_width,fused_scan_record(:,:,1)','edgecolor','none');
    axis equal
    xlim([1 map_width]);ylim([1 map_width]) 
    xlabel('x (m)');ylabel('y (m)');zlabel('m_1(x)');
    view(10,80);
    view(0,90);box on;    
    colormap(jet)
    c=colorbar('eastoutside');
    ylabel(c, 'last visited time (s)')
    title([num2str(cumu_covg),'% of the map explored at t = ',num2str(T-0.2),'s']);        
    
    pause(0.01)
    
    if cumu_covg==100
        break;
    else
        counter = counter+1; 
    end
end