%%
%% LEACH protocol for WSN
%% EE5132/EE5024 CK Tham, ECE NUS
%% Based on code by Smaragdakis, Matta and Bestavros, BU, USA
%%

clear;
IniEng=0.5;
number_of_nodes=100; % Number of Node
num_of_nodes=number_of_nodes;
number_of_rounds=4000; % Number of Round

images_folder = "bcdcp_save_images";

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Field Dimensions - x and y maximum (in meters)
xm_set = [10 30 50 70 90 110 130 150 170 190 200];
ym_set = [10 30 50 70 90 110 130 150 170 190 200];

%x and y coordinates of the BS(will multiply with xm and ym)
%they are relative values, not absolute values
x_BS_location = 1.5;
y_BS_location = 1.5;

%Number of Nodes in the field
n=number_of_nodes;

%Optimal Election Probability of a node
%to become cluster head
p=0.1;

%Energy Model (all values in Joules)
%Initial Energy 
Eo=IniEng;

%Eelec=Etx=Erx
ETX=50*0.000000001;
ERX=50*0.000000001;

%Transmit Amplifier types (Epsilon_{amp})
Efs=10*0.000000000001;
Emp=0.0013*0.000000000001;

do=sqrt(Efs/Emp);

% Data Aggregation Energy
EDA=5*0.000000001;

%maximum number of rounds
rmax=number_of_rounds;

for v = 1:length(xm_set)
    xm = xm_set(v);
    ym = ym_set(v);

    sink.x=x_BS_location*xm_set(v);
    sink.y=y_BS_location*ym_set(v);

rng(0)

%%%%%%%%%%%%%%%%%%%%%%%%% END OF PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%
%Creation of the random Sensor Network
%figure(1);
for i=1:1:n
    S(i).xd=rand(1,1)*xm;
    S(i).yd=rand(1,1)*ym;
    S(i).id=i;
    S(i).G=0; % they were not cluster heads
    %initially there are no cluster heads only nodes
    S(i).type='N'; % nodes
    S(i).E=Eo; % initial energy
    %plot(S(i).xd,S(i).yd,'o');
    %hold on;
end

S(n+1).xd=sink.x;
S(n+1).yd=sink.y;
% plot(S(n+1).xd,S(n+1).yd,'x');
        
%First Iteration
%figure(1);

%counter for CHs
countCHs=0;
%counter for CHs per round
rcountCHs=0;

countCHs;
rcountCHs=rcountCHs+countCHs;
flag_first_dead=0;
PacketRecEnergy     = 0  ;
AggPacketRecvEnergy = [] ;

flag_first_dead=0;
flag_10_percent_dead= 0;
flag_50_percent_dead= 0;
flag_70_percent_dead= 0;
flag_all_dead=0;


round_at_first_dead = 0;
round_at_10_percent_dead = 0;
round_at_50_percent_dead = 0;
round_at_70_percent_dead = 0;
round_at_all_dead = 0;

w_10_percent_dead = 0.1*n;
w_50_percent_dead = 0.5*n;
w_70_percent_dead = 0.7*n;

alive = n;

DM = createDistanceMatrix(S)

%counter for packets transmitted to Base Stations and to Cluster Heads
packets_TO_BS=0;
packets_TO_CH=0;

%%%%%%%%%%%%%%%%%%%%%% START OF MAIN LOOP for rmax rounds %%%%%%%%%%%%%%%%%%
for r=1:1:rmax
  r

hold off;

%Number of dead nodes
dead=0;
alive=0;


%counter for packets transmitted to Base Stations and to Cluster Heads 
%per round
PACKETS_TO_CH(r)=0;
PACKETS_TO_BS(r)=0;

figure(1);

for i=1:1:n
    %checking if there is a dead node
    if (S(i).E<=0)
        plot(S(i).xd,S(i).yd,'red .');
        dead=dead+1;
        
        hold on;    
    end
    if S(i).E>0
        S(i).type='N';
        alive=alive+1;
        plot(S(i).xd,S(i).yd,'o');
        hold on;
    end
end
plot(S(n+1).xd,S(n+1).yd,'x');

STATISTICS.DEAD(r)=dead;
STATISTICS.alive(r)=alive; 
DEAD(r)=dead;


if (dead == n) && (flag_all_dead == 0)
    round_at_all_dead =r;
    flag_all_dead = 1;
end
if (dead >= w_70_percent_dead) && (flag_70_percent_dead == 0)
    round_at_70_percent_dead =r;
    flag_70_percent_dead = 1;
end 
if (dead >= w_50_percent_dead) && (flag_50_percent_dead == 0)
    round_at_50_percent_dead =r;
    flag_50_percent_dead = 1;
end
if (dead >= w_10_percent_dead) && (flag_10_percent_dead == 0)
    round_at_10_percent_dead =r;
    flag_10_percent_dead = 1;
end
if (dead >= 1) && (flag_first_dead == 0)
    round_at_first_dead =r;
    flag_first_dead = 1;
end

countCHs=0;
cluster=1;

% === BCDCP Balanced 2-Clustering ===

% 1) Select nodes that are above average energy
% Get the average energy of all nodes
E_sum = 0;
for i=1:1:n
    if S(i).E > 0
        E_sum = E_sum + S(i).E;
    end
end
E_avg = E_sum / n;
STATISTICS.Total_energy(r)=E_sum;
STATISTICS.Avg_energy(r)=E_avg;

% Note down the nodes that are above average energy
temp_S = [];
for i=1:1:n
    if S(i).E > 0 && S(i).E - E_avg > - E_avg * 0.000000001 % multiply by tolerance
        temp_S = [temp_S S(i)];
    end
end

% Special case for one node only
if length(temp_S) == 1
    % send 1 packet directly to BS
    % Node energy
    node = temp_S(1).id
    distance = sqrt((S(n+1).xd - S(node).xd)^2 + (S(n+1).yd - S(node).yd)^2);

    % might need to change based on energy model
    if distance > do
        S(node).E = S(node).E - ( ETX*2000 + Emp*2000*(distance)^4 );
    else
        S(node).E = S(node).E - ( ETX*2000 + Efs*2000*(distance)^2 );
    end
    packets_TO_BS=packets_TO_BS+1;
    STATISTICS.COUNTCHS(r) = 1;

% Special case for no nodes
elseif length(temp_S) == 0
    STATISTICS.COUNTCHS(r) = 0;
    
% General case for more than one nodes
else
    
    % 2) Balanced Iterative Clustering
    % create array of clusters and cluster heads
    % define nch
    clusters = {temp_S};
    cluster_heads = {randsample(temp_S,1)};
    nch = p * n;
    m_size = n / nch;
    
    % Iterate until the length of clusters is equal to nch
    while length(clusters) < nch
        % select largest cluster and split it to two
        j = 1;
        size = 0;
        for k=1:1:length(clusters)
            if length(clusters{k}) > size
                j = k;
                size = length(clusters{k});
            end
        end

        % exit the loop if largest cluster only contains 1 node
        if size < 2
            m_size = size;
            break
        end

        % Split clusters{j} into two, with:
        % c1 - cluster 1
        % ch1 - cluster head 1
        % c2 - cluster 2
        % ch2 - cluster head 2
        
        [c1, ch1, c2, ch2] = twoClustering(clusters{j}, DM);
        clusters{j} = c1;
        cluster_heads{j} = ch1;
        clusters = [clusters {c2}];
        cluster_heads = [cluster_heads {ch2}];
    end

    STATISTICS.COUNTCHS(r) = length(cluster_heads);



    % Plotting cluster heads as black
    hold on
    for k=1:1:length(cluster_heads)
        plot(S(cluster_heads{k}.id).xd,S(cluster_heads{k}.id).yd,'k*');
        countCHs=countCHs+1;
    end
    STATISTICS.count_CHs(r)=countCHs;
    
    
    % 3) Get shortes route to Base station for each cluster head
    
    % Convert CHs from points cartesian plane into nodes in a graph
    % Convert euclidean distances to weighted edges in a graph
    % Convert point indexes to string, as keys to nodes
    chs = length(cluster_heads);
    s = [];
    t = [];
    w = zeros(chs*(chs + 1)/2);
    edge_count = 1;
    for i=1:1:length(cluster_heads)-1
        for k=i+1:1:length(cluster_heads)
            s = [s string(cluster_heads{i}.id)];
            t = [t string(cluster_heads{k}.id)];
            w(edge_count) = DM(cluster_heads{i}.id, cluster_heads{k}.id);
            edge_count = edge_count+1;
        end
    end
    edge_count = edge_count-1;
    
    % Create a matlab graph object
    G = graph(s(1:edge_count), t(1:edge_count), w(1:edge_count));
    
    % Generate the minimum spanning tree
    [T, pred] = minspantree(G);

    
    % Plot the minimum spanning tree
    %hold on
    for i=1:1:height(T.Edges)
        edge = T.Edges{i,1};
        e_start = str2double(edge{1});
        e_end = str2double(edge{2});

        plot([S(e_start).xd;S(e_end).xd], [S(e_start).yd;S(e_end).yd]);

    end

    % ========= CH-to-CH route & energy calculation =========

    % random pick of main cluster head sending to Base Station
    main_ch = randsample(cluster_heads,1);
    main_ch = main_ch{1}.id;

    text(S(main_ch).xd,S(main_ch).yd,'M', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
    plot([S(main_ch).xd;S(n+1).xd], [S(main_ch).yd;S(n+1).yd]);

    % Iterate for each cluster head
    for j=1:1:length(cluster_heads)
        
        % ===== receive from cluster nodes =====
        
        % Get cluster head id
        ch = cluster_heads{j}.id;

        % Iterate through each nodes in the cluster
        for k=1:1:length(clusters{j})

            % If nodes are not dead
            if clusters{j}(k).E > 0 

                % Transfer from node to ch
                if clusters{j}(k).id ~= ch
                    
                    % Node id and energy
                    node = clusters{j}(k).id;
                    distance = DM(ch, node);

                    % Transmission energy dissipation
                    if distance > do
                        S(node).E = S(node).E - ( ETX*2000 + Emp*2000*(distance)^4 );
                    else
                        S(node).E = S(node).E - ( ETX*2000 + Efs*2000*(distance)^2 );
                    end

                    % Receive energy dissipation
                    S(ch).E = S(ch).E - (ERX + EDA) * 2000;
                else
                    1;
                end
            end
        end

        % ===== ch-to-ch routing =====
        
        % Find the shortest path in the minimum spanning tree
        ch_path = shortestpath(T, string(ch), string(main_ch));


        % For every node in the cluster path
        for k=1:1:length(ch_path)-1
            
            % Get current node and the next node
            node = str2double(ch_path{k});
            node_next = str2double(ch_path{k+1});

            
            % Transmission energy dissipation
            distance = DM(node, node_next);
            
            S(node).E = S(node).E - ( ETX*2000 + Emp*2000*(distance*distance) );

            % Receive energy dissipation
            S(node_next).E = S(node_next).E - (ERX + EDA) * 2000;
        end

        % main_ch send data packet to Base Station
        distance = DM(n+1, main_ch);

        % Transmission energy dissipation
        if distance > do
            S(node).E = S(node).E - ( ETX*2000 + Emp*2000*(distance)^4 );
        else
            S(node).E = S(node).E - ( ETX*2000 + Efs*2000*(distance)^2 );
        end

        packets_TO_BS=packets_TO_BS+1;

    end

end

packets_TO_BS;
PACKETS_TO_BS(r) = packets_TO_BS;
STATISTICS.PACKETS_TO_BS(r)=packets_TO_BS;

countCHs;
rcountCHs=rcountCHs+countCHs;

PacketRecEnergy = 0 ;

for i=1:n
    if(S(i).E < 0)
        temp = 0 ;
    else
        temp = S(i).E ;
    end
    PacketRecEnergy = PacketRecEnergy + Eo - temp; % dissipated energy

end

AggPacketRecvEnergy = [AggPacketRecvEnergy PacketRecEnergy] ;

end
%%%%%%%%%%%%%%%%%%%%%% END OF MAIN LOOP for rmax rounds %%%%%%%%%%%%%%%%%%

energycontainer(v)=AggPacketRecvEnergy(rmax);
STATISTICS2=STATISTICS;
FD2=round_at_first_dead;
TD2=round_at_10_percent_dead;
AD2=round_at_all_dead;

% plot cluster head num of every round
r=0:rmax;
% -------------------containers for plot----------------------------------
container_for_plot(v).first = STATISTICS2.PACKETS_TO_BS*10;
container_for_plot(v).second = AggPacketRecvEnergy;
container_for_plot(v).third = STATISTICS2.alive;
container_for_plot(v).fourth = STATISTICS2.COUNTCHS;
container_for_plot(v).dead_round = [round_at_first_dead, round_at_10_percent_dead,round_at_50_percent_dead, round_at_70_percent_dead,round_at_all_dead];

container_for_plot(v).lastIdx = max(find(container_for_plot(v).third == num_of_nodes));


end

% ++++++++++++++++++++ Plotting Figure 4(f) Different network size
minIdx = 2001;
for i=1:length(container_for_plot)
    minIdx = min(container_for_plot(i).lastIdx, minIdx);
end
minIdx

energycontainer = [];
for i=1:length(container_for_plot)
    array = container_for_plot(i).second;
    energycontainer = [energycontainer array(minIdx)];
end

maximum = max(energycontainer);
normalised = energycontainer/ maximum;

fig =figure;
% plot(p_set*100,energycontainer)
plot(xm_set,normalised,'LineWidth',2);
title_label = 'Different Network Size vs Normalised Energy Dissipation';
title([title_label])

xlabel('Different NetWork Diameter')
ylabel('Normalised Energy Dissipation')
file_zero = images_folder + "/Different Network Size vs Normalised Energy Dissipation.png";
saveas(fig,[file_zero]);

% -------------------------Plotting all combined--------------------------

% different types of lines for plotting
line_types = char('-r','-g','-b','-c','-m','-y','-k','--r','--g','--b','--c','--m','--y','--k');
% legend
legend_f1 = 'D-size = 10';
legend_f2 = 'D-size = 30';
legend_f3 = 'D-size = 50';
legend_f4 = 'D-size = 70';
legend_f5 = 'D-size = 90';
legend_f6 = 'D-size = 110';
legend_f7 = 'D-size = 130';
legend_f8 = 'D-size = 150';
legend_f9 = 'D-size = 170';
legend_f10 = 'D-size = 190';
legend_f11 = 'D-size = 200'; 

%--------------------------------------------------------------
first = figure;
for i=1:length(container_for_plot)
    plot(container_for_plot(i).first,line_types(i,:),'LineWidth',1.5);
    hold on
end
legend(legend_f1,legend_f2,legend_f3,legend_f4,legend_f5,legend_f6,legend_f7,legend_f8,legend_f9,legend_f10,legend_f11,'Location','northwest');
title_label = 'Plot for Different Network Size values';
title([title_label])

xlabel('Time(rounds) (Fig 7(a) TWC)')
ylabel('Number of data signals received at the base station')
file_time = images_folder + "/data pkts at BS vs. rounds.png";
saveas(first,[file_time]);
%--------------------------------------------------------------
second = figure;
for i=1:length(container_for_plot)
    plot(container_for_plot(i).second,container_for_plot(i).first,line_types(i,:),'LineWidth',1.5);
    hold on
end
legend(legend_f1,legend_f2,legend_f3,legend_f4,legend_f5,legend_f6,legend_f7,legend_f8,legend_f9,legend_f10,legend_f11,'Location','northwest');
title_label = 'Plot for Different Network Size values';
title([title_label])

xlabel('Energy(J)(Fig 7(b) TWC)')
ylabel('Number of data signals received at the base station')
file_energy = images_folder + "/data pkts at BS vs. energy.png";
saveas(second,[file_energy]);
%--------------------------------------------------------------
Third = figure;
for i=1:length(container_for_plot)
    plot(container_for_plot(i).third,line_types(i,:),'LineWidth',1.5);
    hold on
end
legend(legend_f1,legend_f2,legend_f3,legend_f4,legend_f5,legend_f6,legend_f7,legend_f8,legend_f9,legend_f10,legend_f11,'Location','southwest');
title_label = 'Plot for Different Network Size values';
title([title_label])

xlabel('Time(rounds) (Fig 8(a) TWC)')
ylabel('Number of nodes alive')
file_time_nodes_alive = images_folder + "/alive nodes vs. rounds.png";
saveas(Third,[file_time_nodes_alive]);

%--------------------------------------------------------------
fourth = figure;
for i=1:length(container_for_plot)
    plot(container_for_plot(i).first,container_for_plot(i).third,line_types(i,:),'LineWidth',1.5);
    hold on
end
legend(legend_f1,legend_f2,legend_f3,legend_f4,legend_f5,legend_f6,legend_f7,legend_f8,legend_f9,legend_f10,legend_f11,'Location','southwest');
title_label = 'Plot for Different Network Size values';
title([title_label])

xlabel('Number of data items received at BS (Fig 8(b) TWC)')
ylabel('Number of nodes alive')

file_dataitem_nodes_alive = images_folder + "/alive nodes vs. data pkts at BS.png";
saveas(fourth,[file_dataitem_nodes_alive]);

%--------------------------------------------------------------
fifth = figure;
for i=1:length(container_for_plot)
    plot(cumsum(container_for_plot(i).fourth),line_types(i,:),'LineWidth',1.5);
    hold on
end
legend(legend_f1,legend_f2,legend_f3,legend_f4,legend_f5,legend_f6,legend_f7,legend_f8,legend_f9,legend_f10,legend_f11,'Location','northwest');
title_label = 'Plot for Different Network Size values COUNTCH';
title([title_label])

xlabel('Time (rounds)')
ylabel('Number of accumulative cluster heads')

file_countch = images_folder + "/accumulative CH vs. data pkts at BS.png";
saveas(fifth,[file_countch]);

%---------------------------------bar figure--------------------------------
% x axis: different percentage of death. each category(each percentage)'s bars
% correspond to different diameters
first_death_for_diff_xm = [];
ten_death_for_diff_xm = [];
fifty_death_for_diff_xm = [];
seventy_death_for_diff_xm = [];
all_death_for_diff_xm = [];
for i=1:length(container_for_plot)
    first_death_for_diff_xm = [first_death_for_diff_xm container_for_plot(i).dead_round(1)];
    ten_death_for_diff_xm = [ten_death_for_diff_xm container_for_plot(i).dead_round(2)];
    fifty_death_for_diff_xm = [fifty_death_for_diff_xm container_for_plot(i).dead_round(3)];
    seventy_death_for_diff_xm = [seventy_death_for_diff_xm container_for_plot(i).dead_round(4)];
    all_death_for_diff_xm = [all_death_for_diff_xm container_for_plot(i).dead_round(5)];   
end
% names are too long, rename the variables
d1 = first_death_for_diff_xm;
d10 = ten_death_for_diff_xm;
d50 = fifty_death_for_diff_xm;
d70 = seventy_death_for_diff_xm;
d100 = all_death_for_diff_xm;

barfig = figure;
x = categorical({'1st','10%','50%','70%','100%'});
x = reordercats(x,{'1st','10%','50%','70%','100%'});
y=[d1(1) d1(2) d1(3) d1(4) d1(5) d1(6) d1(7) d1(8) d1(9) d1(10) d1(11);d10(1) d10(2) d10(3) d10(4) d10(5) d10(6) d10(7) d10(8) d10(9) d10(10) d10(11);d50(1) d50(2) d50(3) d50(4) d50(5) d50(6) d50(7) d50(8) d50(9) d50(10) d50(11);d70(1) d70(2) d70(3) d70(4) d70(5) d70(6) d70(7) d70(8) d70(9) d70(10) d70(11);d100(1) d100(2) d100(3) d100(4) d100(5) d100(6) d100(7) d100(8) d100(9) d100(10) d100(11)]; 
bar(x,y,'group');
title_label='The round when specific percentage of nodes die for different diameters';
title([title_label]);
xlabel('Amount of deaths');
ylabel('Number of rounds');

figure_name = [images_folder + "/death round.png"];
saveas(barfig,[figure_name]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   STATISTICS    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                     %
%  DEAD  : a rmax x 1 array of number of dead nodes/round 
%  CLUSTERHS : a rmax x 1 array of number of Cluster Heads/round
%  PACKETS_TO_BS : a rmax x 1 array of number packets send to Base Station/round
%  PACKETS_TO_CH : a rmax x 1 array of number of packets send to ClusterHeads/round
%  first_dead: the round where the first node died                   
%                                                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function DM = createDistanceMatrix(S)
    DM = zeros(length(S), length(S))
    for i=1:1:length(S)-1
        for k=i:1:length(S)
            dist_temp = sqrt((S(i).xd - S(k).xd)^2 + (S(i).yd - S(k).yd)^2);
            DM(i,k) = dist_temp;
            DM(k,i) = dist_temp;
        end
    end
end


function [c1,ch1,c2,ch2] = twoClustering(c, DM)

    dist = -1;
    ch1 = -1;
    ch2 = -1;
    % find largest distance
    for i=1:1:length(c)-1
        for k=i+1:1:length(c)
            
            dist_temp = DM(i,k);
            if dist_temp > dist
                dist = dist_temp;
                ch1 = c(i);
                ch2 = c(k);
            end
        end
    end
    
    c1 = [];
    c2 = [];
    
    % group nodes to the two clusters
    for i=1:1:length(c)
        d_ch1 = DM(c(i).id, ch1.id);
        d_ch2 = DM(c(i).id, ch2.id);
        
        if d_ch1 < d_ch2
            c1 = [c(i) c1];
            
        else
            c2 = [c(i) c2];
            
        end
    end
    
    % fix imbalances
    diff = length(c1) - length(c2);
    if diff ~= 0
        
        % if |c1| > |c2| find nodes in c1 closest to c2 then move them
        if diff > 0
            tomove = floor(diff / 2);
            dists = [];
            for i=1:1:length(c1)
                
                dists = [dists, DM(c1(i).id, ch2.id)];
            end
            dists;
            [X,idx] = sort(dists);
            
            c2 = [c1(idx(1:tomove)) c2];
            c1 = c1(idx(tomove+1:length(idx)));
        
        % if |c1| < |c2| find nodes in c2 closest to c1 then move them
        elseif diff < 0
            tomove = - floor(diff / 2);
            dists = [];
            for i=1:1:length(c2)
                
                dists = [dists, DM(c2(i).id, ch1.id)];
            end
            dists;
            [X,idx] = sort(dists);
            
            c1 = [c2(idx(1:tomove)) c1];
            c2 = c2(idx(tomove+1:length(idx)));
        end
    end

end

