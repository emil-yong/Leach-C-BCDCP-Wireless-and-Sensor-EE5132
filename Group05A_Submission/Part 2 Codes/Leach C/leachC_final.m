clear
clc

%%
%% LEACH C protocol for WSN
%% Done by Group 5 EE5132 NUS CA2 Assignment

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameters used is similar to adaptive Leach provided

IniEng=0.5;%0.5; % Initial Energy of Every Node
num_of_nodes=100; % Number of Node
number_of_rounds=4000; % Number of Round


pkt_size = 2000;

%Field Dimensions - x and y maximum (in meters)
% Our group chose to vary the size of network diameter
xm_set = [10 30 50 70 90 110 130 150 170 190 200];
ym_set = [10 30 50 70 90 110 130 150 170 190 200];


%x and y coordinates of the BS(will multiply with xm and ym)
% Relative values to determined the location of BS 
% For differnt network size
x_BS_location = 1.5;
y_BS_location = 1.5;

%Number of Nodes in the field
n=num_of_nodes;

p=0.1;


%Energy Model (all values in Joules)
%Initial Energy 
Eo=IniEng;


% Efs is the epsilon for Free Space
% Emp is the Epsilon for Multipath Fading
% Note that Experiment set up by Heinzelman 2002 uses a simple model

%Eelec=Etx=Erx
ETX=50*0.000000001;
ERX=50*0.000000001;
%Transmit Amplifier types
Efs=10*0.000000000001;
Emp=0.0013*0.000000000001;

%Data Aggregation Energy
EDA=5*0.000000001;

varying_Eo=0;

%maximum number of rounds
rmax=number_of_rounds;

do=sqrt(Efs/Emp);

% Initialize to keep track of the max and min to the BS
% to find the K_optimal

for v = 1:length(xm_set)
    xm = xm_set(v);
    ym = ym_set(v);

    sink.x=x_BS_location*xm_set(v);
    sink.y=y_BS_location*ym_set(v);

    d_to_BS_min = 10000000000;
    d_to_BS_max = 0;
    d_to_BS_avg = 0;

for i=1:1:n
    S(i).xd=rand(1,1)*xm;
    S(i).yd=rand(1,1)*ym;
    S(i).G=0;
    % Setting the initial energy to be Eo
    % can play around with a but we set it to be 0
    % for consistancy with the leach
    S(i).E=Eo*(1+rand*varying_Eo);
    %initially there are no cluster heads only nodes
    S(i).type='N';
    %compute distance between this node and BS
    distance=sqrt( (S(i).xd-(sink.x) )^2 + (S(i).yd-(sink.y) )^2 );
    % keep track on the distance to all the nodes to BS
    % Find the eventual d_to_bs_min
    % Find the eventual d_to_bs_max
    d_to_BS_avg = d_to_BS_avg + distance;
    d_to_BS_min = min(d_to_BS_min,distance);
    d_to_BS_max = max(d_to_BS_max,distance);
end


% -----compute optimal k (k_opt)-----
N = num_of_nodes;
M = xm;

k_opt_avg = sqrt(N/(2*pi))*do*(M/((d_to_BS_avg/N)^2));
k_opt_min = sqrt(N/(2*pi))*do*(M/(d_to_BS_max^2));
k_opt_max = sqrt(N/(2*pi))*do*(M/(d_to_BS_min^2));

if k_opt_max > N
    k_opt_max = N;
end

xm
d_to_BS_min 
d_to_BS_max 
d_to_BS_avg 
k_opt_avg 
k_opt_min
k_opt_max 
% -----compute optimal k (k_opt)-----

container_for_plot(v).koptmin = k_opt_min;
container_for_plot(v).koptmax = k_opt_max;

%BS store in variable S
S(n+1).xd=sink.x;
S(n+1).yd=sink.y;

%Intialised count for CHs
countCHs=0;
cluster=1;

%Keep track on different percentage of death
flag_first_dead=0;
flag_10_percent_dead= 0;
flag_50_percent_dead= 0;
flag_70_percent_dead= 0;
flag_all_dead=0;

% Find the rounds at which the death occurs for each percentage
round_at_first_dead = 0;
round_at_10_percent_dead = 0;
round_at_50_percent_dead = 0;
round_at_70_percent_dead = 0;
round_at_all_dead = 0;

% different percentage of death
w_10_percent_dead = 0.1*n;
w_50_percent_dead = 0.5*n;
w_70_percent_dead = 0.7*n;

alive=n;
%counter for bit transmitted to Bases Station and to Cluster Heads
packets_TO_BS=0;
packets_TO_CH=0;

PacketRecEnergy     = 0  ;
AggPacketRecvEnergy = [] ;


for r=0:1:rmax     
    %r
    if(mod(r, round(1/p))==0)
        for i=1:1:n
            S(i).G=0; 
            %S(i).cl=0;
        end
    end

    dead=0;
    for i=1:1:n
        if (S(i).E<=0)
            dead=dead+1;
            % switch case to keep track on the which round
            % percentage of death occurs
            switch dead
                case 1
                    if (flag_first_dead == 0)
                        round_at_first_dead =r;
                        flag_first_dead = 1;
                    end
            
                case w_10_percent_dead
                    if (flag_10_percent_dead == 0)
                        round_at_10_percent_dead =r;
                        flag_10_percent_dead = 1;
                    end     
                case w_50_percent_dead
                    if (flag_50_percent_dead == 0)
                        round_at_50_percent_dead =r;
                        flag_50_percent_dead = 1;
                    end                            
                case w_70_percent_dead
                    if (flag_70_percent_dead == 0)
                        round_at_70_percent_dead =r;
                        flag_70_percent_dead = 1;
                    end 
                case n
                    if (flag_all_dead == 0)
                        round_at_all_dead =r;
                        flag_all_dead = 1;
                    end 
            end
        end
        if S(i).E>0
            S(i).type='N';
        end
    end
    STATISTICS.DEAD(r+1)=dead;
    STATISTICS.alive(r+1)=alive-dead;  

    %Calculate the average energy of the whole network%%%%%%%%%%%%%%%%%%%%
    Etotal=0;
    for i=1:n
        if S(i).E>0
            Etotal=Etotal+S(i).E;
        end
    end
    Eavg=Etotal/n;
    STATISTICS.Total_energy(r+1)=Etotal;
    %% Total energy remaining
    STATISTICS.AvgEnergy(r+1)=Eavg;

    countCHs=0;
    cluster=1;
    for i=1:1:n
        if(S(i).E>0)
            temp_rand=rand;
            if ( (S(i).G)<=0)  
                if (temp_rand<= (p/(1-p*mod(r,round(1/p)))))
                    if (S(i).E>Eavg) % should be a CH
                        countCHs=countCHs+1;
                        packets_TO_BS=packets_TO_BS+1;
                        PACKETS_TO_BS(r+1)=packets_TO_BS; %independent variable
                        S(i).type='C';
                        S(i).G=round(1/p)-1;
                        C(cluster).xd=S(i).xd;
                        C(cluster).yd=S(i).yd;
                        distance=sqrt( (S(i).xd-(S(n+1).xd) )^2 + (S(i).yd-(S(n+1).yd) )^2 );
                        C(cluster).distance=distance;
                        C(cluster).id=i;
                        X(cluster)=S(i).xd;
                        Y(cluster)=S(i).yd;
                        cluster=cluster+1;
                        
                        distance;
                        if (distance>do) %！！！！equation 9
                            S(i).E=S(i).E- ( (ETX+EDA)*(pkt_size) + Emp*pkt_size*( distance*distance*distance*distance )); 
                        end
                        if (distance<=do)
                            S(i).E=S(i).E- ( (ETX+EDA)*(pkt_size)  + Efs*pkt_size*( distance * distance )); 
                        end
                    end
                end     
            end
              
        end 
    end
    STATISTICS.COUNTCHS(r+1)=countCHs;
    %pause;

    for i=1:1:n
        if ( S(i).type=='N' && S(i).E>0 ) % if it is a node OR alive
            if(cluster-1>=1)
                min_dis=sqrt( (S(i).xd-S(n+1).xd)^2 + (S(i).yd-S(n+1).yd)^2 );
                min_dis_cluster=1;
                for c=1:1:cluster-1
                    temp=min(min_dis,sqrt( (S(i).xd-C(c).xd)^2 + (S(i).yd-C(c).yd)^2 ) );
                    if ( temp<min_dis )
                        min_dis=temp;
                        min_dis_cluster=c;
                    end
                end
                 
                if(min_dis_cluster~=1)  % 
                    min_dis;
                    if (min_dis>do)
                        S(i).E=S(i).E- ( ETX*(pkt_size) + Emp*pkt_size*( min_dis * min_dis * min_dis * min_dis)); 
                    end
                    if (min_dis<=do)
                        S(i).E=S(i).E- ( ETX*(pkt_size) + Efs*pkt_size*( min_dis * min_dis)); 
                    end

                    % Because CH is the intermediate node to the BS
                        % Data received from other nodes also will have
                        % energy dissapation.
                        % Data Aggregation etc
                    S(C(min_dis_cluster).id).E = S(C(min_dis_cluster).id).E- ( (ERX + EDA)*pkt_size ); 
                    packets_TO_CH=packets_TO_CH+1; 
                    % dont need to consider the amount of packet send
                        % to the BS in this for loop because when selecting
                        % the cluster head, it is already considered.
                    

                else % Closer to the BS than all other cluster head. Send to BS straight
                    min_dis;
                    if (min_dis>do)
                        S(i).E=S(i).E- ( ETX*(pkt_size) + Emp*pkt_size*( min_dis * min_dis * min_dis * min_dis)); 
                    end
                    if (min_dis<=do)
                        S(i).E=S(i).E- ( ETX*(pkt_size) + Efs*pkt_size*( min_dis * min_dis)); 
                    end
                    packets_TO_BS=packets_TO_BS+1;    % send pkt to BS
                    
                end
                S(i).min_dis=min_dis;
                S(i).min_dis_cluster=min_dis_cluster;
            else 
                    % cluster-1<1 consider where there is only 1 cluster
                    % in the code provided by Leach.m, the code only
                    % include where there is always more than 1 cluster
                    % if there only exisit one cluster means base
                    % station is considered as the clusterhead for all nodes to transmit 
                min_dis=sqrt( (S(i).xd-S(n+1).xd)^2 + (S(i).yd-S(n+1).yd)^2 );
                if (min_dis>do)
                    S(i).E=S(i).E- ( ETX*(pkt_size) + Emp*pkt_size*( min_dis * min_dis * min_dis * min_dis)); 
                end
                if (min_dis<=do)
                    S(i).E=S(i).E- ( ETX*(pkt_size) + Efs*pkt_size*( min_dis * min_dis)); 
                end
                packets_TO_BS=packets_TO_BS+1;
            end
        end
    end
    STATISTICS.PACKETS_TO_CH(r+1)=packets_TO_CH;
    STATISTICS.PACKETS_TO_BS(r+1)=packets_TO_BS;

    PacketRecEnergy = 0 ;
    
    for i=1:n
        if(S(i).E < 0)
            temp = 0 ;
        else
            temp = S(i).E ;
        end
            PacketRecEnergy = PacketRecEnergy + Eo - temp; 
            % total energy used accumulated
    end
    AggPacketRecvEnergy = [AggPacketRecvEnergy PacketRecEnergy] ;
end

energycontainer(v)=AggPacketRecvEnergy(number_of_rounds);

%energyplt(v)=AggPacketRecvEnergy;

%4 return values:
STATISTICS_2_store=STATISTICS;
% plot cluster head num of every round
r=0:number_of_rounds;
% -------------------containers for plot----------------------------------
container_for_plot(v).first = STATISTICS_2_store.PACKETS_TO_BS;
container_for_plot(v).second = AggPacketRecvEnergy;
container_for_plot(v).third = STATISTICS_2_store.alive;
container_for_plot(v).fourth = STATISTICS_2_store.COUNTCHS;
container_for_plot(v).dead_round = [round_at_first_dead, round_at_10_percent_dead,round_at_50_percent_dead, round_at_70_percent_dead,round_at_all_dead];


container_for_plot(v).lastIdx = max(find(container_for_plot(v).third == num_of_nodes));


% ----------------------------KOPT ----------------------------
fig = figure;
plot(r,STATISTICS.COUNTCHS,'y',r,repelem(k_opt_min,number_of_rounds+1),'r',r,repelem(k_opt_max,number_of_rounds+1),'b', ...
    r,repelem(k_opt_avg,number_of_rounds+1),'g','LineWidth',1);
legend('CH num','min kopt','max kopt','avg kopt','Location','northeast');

xlabel('round');
ylabel('No of cluster heads');
title(['No of cluster heads over time for diameter=' num2str(xm)]);

saveas(fig,['leachC_save_images/' num2str(xm) 'kopt.png']);
end



% ++++++++++++++++++++ Plotting Figure 4(f) Different network size
minIdx = number_of_rounds+1;
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
plot(xm_set,normalised,'LineWidth',2);
title_label = 'Different Network Size vs Normalised Energy Dissipation';
title([title_label])

xlabel('Different NetWork Diameter')
ylabel('Normalised Energy Dissipation')
file_zero = 'leachC_save_images/Different Network Size vs Normalised Energy Dissipation.png';
saveas(fig,[file_zero]);

% -------------------------Plotting all combined--------------------------

% different types of lines for plotting
line_types = char('-r','-g','-b','-c','-m','-y','-k','--r','--g','--b','--c','--m','--y','--k');
% legend
legend_f1 = 'D-size= 10';
legend_f2 = 'D-size= 30';
legend_f3 = 'D-size =50';
legend_f4 = 'D-size =70';
legend_f5 = 'D-size =90';
legend_f6 = 'D-size= 110';
legend_f7 = 'D-size= 130';
legend_f8 = 'D-size =150';
legend_f9 = 'D-size =170';
legend_f10 = 'D-size =190';
legend_f11 = 'D-size =200'; 

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
file_time = 'leachC_save_images/data pkts at BS vs. rounds.png';
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
file_energy = 'leachC_save_images/data pkts at BS vs. energy.png';
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
file_time_nodes_alive = 'leachC_save_images/alive nodes vs. rounds.png';
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

file_dataitem_nodes_alive = 'leachC_save_images/alive nodes vs. data pkts at BS.png';
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

file_countch = 'leachC_save_images/accumulative CH vs. data pkts at BS.png';
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

figure_name = ['leachC_save_images/death round.png'];
saveas(barfig,[figure_name]);
