%%
%% LEACH protocol for WSN
%% EE5132/EE5024 CK Tham, ECE NUS
%% Based on code by Smaragdakis, Matta and Bestavros, BU, USA
%%

clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Field Dimensions - x and y maximum (in meters)
% Our group chose to vary the size of network diameter
xm_set = [10 30 50 70 90 110 130 150 170 190 200];
ym_set = [10 30 50 70 90 110 130 150 170 190 200];


%Number of Nodes in the field
n=100;

%Optimal Election Probability of a node
%to become cluster head
p=0.1;


%Energy Model (all values in Joules)
%Initial Energy 
Eo=0.5;

% Efs is the epsilon for Free Space
% Emp is the Epsilon for Multipath Fading
% Note that Experiment set up by Heinzelman 2002 uses a simple model

%Eelec=Etx=Erx
ETX=50*0.000000001;
ERX=50*0.000000001;
%Transmit Amplifier types
Efs=10*0.000000000001;
%Transmit Amplifier types (Epsilon_{amp})
Emp=0.0013*0.000000000001;%Align with Leach C and BCDCP paper
% Data Aggregation Energy
EDA=5*0.000000001;%Align with Leach C and BCDCP paper

do=sqrt(Efs/Emp);

pkt_size = 2000;

%maximum number of rounds
rmax=4000;

for v = 1:length(xm_set)
     xm = xm_set(v);
     ym = ym_set(v);

     %x and y Coordinates of the BS
     % Relative values to determined the location of BS 
     % For differnt network size
     sink.x=1.5*xm_set(v);
     sink.y=1.5*ym_set(v);

%%%%%%%%%%%%%%%%%%%%%%%%% END OF PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%
%Creation of the random Sensor Network
%figure(1);
for i=1:1:n
    S(i).xd=rand(1,1)*xm; % X coordinate of the node
    S(i).yd=rand(1,1)*ym; % Y coordinate of the node 
    S(i).G=0; % set the G to be 0
    %initially there are no cluster heads only nodes
    S(i).type='N';
    S(i).E=Eo;   % E is the energy
    %plot(S(i).xd,S(i).yd,'o');
    %hold on;
end
% ther very last of S varaible is the SINK node which is the information of
% the BASE Station 
S(n+1).xd=sink.x;
S(n+1).yd=sink.y;
%plot(S(n+1).xd,S(n+1).yd,'x');
        
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
r_label = [];

% record round number of specific death percentage
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

%%%%%%%%%%%%%%%%%%%%%% START OF MAIN LOOP for rmax rounds %%%%%%%%%%%%%%%%%%
for r=1:1:rmax %% Create a for loop to iterate all network diameter 
  r % This is to print out the iteration of the loop

  % AFter 1/p round it will reset the nodes that are not cluster head
  % to G  = 0
  %Operation for epoch
  if(mod(r, round(1/p) )==0)
    for i=1:1:n
        S(i).G=0;
    end
  end

%hold off;

%Number of dead nodes
dead=0;

%counter for packets transmitted to Base Stations and to Cluster Heads
packets_TO_BS=0;
packets_TO_CH=0;

%counter for packets transmitted to Base Stations and to Cluster Heads 
%per round
PACKETS_TO_CH(r)=0;
PACKETS_TO_BS(r)=0;

%figure(1);

for i=1:1:n
    %checking if there is a dead node
    if (S(i).E<=0)
        dead=dead+1;
        % if there is a dead node, add the count
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

STATISTICS(r).DEAD=dead;  % Count how many in the rounds are dead
DEAD(r)=dead;

%When the first node dies
if (dead==1)
    if(flag_first_dead==0)
        first_dead=r ;   % record the round where the first dead occurs
        flag_first_dead=1;
    end
end

countCHs=0;
cluster=1;
for i=1:1:n
   if(S(i).E>0) % still alive
     temp_rand=rand;     
     if ((S(i).G)<=0) % is it in set G?

     %Election of Cluster Heads
       if(temp_rand<=(p/(1-p*mod(r,round(1/p)))))
            countCHs=countCHs+1;
            packets_TO_BS=packets_TO_BS+1;
            PACKETS_TO_BS(r)=packets_TO_BS;
            
            S(i).type='C';  % set that node i to to be the CH
            S(i).G=round(1/p)-1; % big value
            C(cluster).xd=S(i).xd;
            C(cluster).yd=S(i).yd;
            %plot(S(i).xd,S(i).yd,'k*'); % plot inside the o with a *
            
            %distance from the cluster head to the base station
            distance=sqrt( (S(i).xd-(S(n+1).xd))^2 + (S(i).yd-(S(n+1).yd))^2 );
            C(cluster).distance=distance;
            C(cluster).id=i;
            cluster=cluster+1;
            
            % Calculation of Energy dissipated
            distance;
            if (distance>do)
                S(i).E=S(i).E- ( (ETX)*(pkt_size) + Emp*pkt_size*( distance*distance*distance*distance )); 
            end
            if (distance<=do)
                S(i).E=S(i).E- ( (ETX)*(pkt_size)  + Efs*pkt_size*( distance * distance )); 
            end
            % remaining of the energy after it minus the 
            % transmitted energy

       end
     end
   end
end

STATISTICS(r).CLUSTERHEADS=cluster-1; %maybe because C(cluster) cannot index 0
CLUSTERHS(r)=cluster-1;

% Election of Associated Cluster Head for Normal Nodes
for i=1:1:n
   if ( S(i).type=='N' && S(i).E>0 )
     if(cluster-1>=1) % iterating through all the cluster head
         %if it is close to BS then u can send to the cluster head 1 already
         % May form a bad cluster
       min_dis=sqrt( (S(i).xd-S(n+1).xd)^2 + (S(i).yd-S(n+1).yd)^2 );
       min_dis_cluster=1;

       for c=1:1:cluster-1
           % if the node is closer to the BS
           % it is assigned to the first cluster
           % Else it is assigned to the shortest CH
           temp=min(min_dis,sqrt( (S(i).xd-C(c).xd)^2 + (S(i).yd-C(c).yd)^2 ) );
           if ( temp<min_dis )
               min_dis=temp;
               min_dis_cluster=c;
           end
       end
       
       % Energy dissipated to associated Cluster Head
       min_dis;
       if (min_dis>do)
           S(i).E=S(i).E- ( ETX*(pkt_size) + Emp*pkt_size*( min_dis * min_dis * min_dis * min_dis)); 
       end
       if (min_dis<=do)
           S(i).E=S(i).E- ( ETX*(pkt_size) + Efs*pkt_size*( min_dis * min_dis)); 
       end
       
       % Energy dissipated at associated Cluster Head
       if(min_dis>0)
           % data received from other nodes
           % there is data aggregation so there is processing

           S(C(min_dis_cluster).id).E = S(C(min_dis_cluster).id).E - (ERX+EDA)*2000; 
           PACKETS_TO_CH(r)=n-dead-cluster+1; 
           % all the packet send to their resepective cluster head
       end
       S(i).min_dis=min_dis;
       S(i).min_dis_cluster=min_dis_cluster;
     end
  end
end
%hold on;

countCHs;
rcountCHs=rcountCHs+countCHs;

PacketRecEnergy = 0 ;

for i=1:n
    if(S(i).E < 0)
        temp = 0 ;
    else
        temp = S(i).E ; % because S(i).E is already the accumulated 
    end
    PacketRecEnergy = PacketRecEnergy + Eo - temp;
    % Accumulated consumed energy used per round from all nodes

end

%append to roundmax 
%roundmax array size
AggPacketRecvEnergy = [AggPacketRecvEnergy PacketRecEnergy] ;
% for every round expenses energy
r_label = [r_label r];


end

energycontainer(v)=AggPacketRecvEnergy(rmax);%/(n*Eo);
energyplt(v).first=AggPacketRecvEnergy;


%%%%%%%%%%%%%%%%%%%%%% END OF MAIN LOOP for rmax rounds %%%%%%%%%%%%%%%%%%
container_for_plot(v).first = cumsum(PACKETS_TO_BS);
container_for_plot(v).second = AggPacketRecvEnergy;
container_for_plot(v).third = (n*ones(1,rmax)-DEAD);%alive
container_for_plot(v).dead_round = [round_at_first_dead, round_at_10_percent_dead,round_at_50_percent_dead, round_at_70_percent_dead,round_at_all_dead];
container_for_plot(v).clusterheads = CLUSTERHS;


container_for_plot(v).lastIdx = max(find(container_for_plot(v).third == n));

% plot figures for each diameter value:
xm_str = num2str(xm);
% % ==============================bar figure===================================
% zero = figure;
% x = categorical({'1st','10%','50%','70%','100%'});
% x = reordercats(x,{'1st','10%','50%','70%','100%'});
% y=[round_at_first_dead round_at_10_percent_dead round_at_50_percent_dead round_at_70_percent_dead round_at_all_dead]; 
% bar(x,y);
% title_label='The round when specific percentage of nodes die, for diameter=';
% title([title_label xm_str]);
% xlabel('Amount of death');
% ylabel('Number of rounds');
% 
% figure_name = ['leach_save_images/',xm_str,' death round.png'];
% saveas(zero,[figure_name]);
% % =================================================================
% first = figure; 
% 
% plot(cumsum(PACKETS_TO_BS))
% 
% title_label = 'Plot for Network Diameter = ';
% title([title_label xm_str])
% xlabel('Time(rounds) (Fig 7(a) TWC)')
% ylabel('Number of data signals received at the base station')
% 
% figure_name = ['leach_save_images/',xm_str,' data items at BS vs. rounds.png'];
% saveas(first,[figure_name]);
% 
% % cumsum is cumulative sum
% % =================================================================
% 
% second = figure;
% plot(AggPacketRecvEnergy,cumsum(PACKETS_TO_BS))
% 
% title_label = 'Plot for Network Diameter = ';
% title([title_label xm_str])
% 
% xlabel('Energy(J)(Fig 7(b) TWC)')
% ylabel('Number of data signals received at the base station')
% 
% figure_name = ['leach_save_images/',xm_str,' data items at BS vs. energy.png'];
% saveas(second,[figure_name]);
% 
% % % =================================================================
% 
% third= figure;
% plot(100*ones(1,rmax)-DEAD)
% 
% title_label = 'Plot for Network Diameter = ';
% title([title_label xm_str])
% 
% xlabel('Time(rounds) (Fig 8(a) TWC)')
% ylabel('Number of nodes alive')
% 
% figure_name = ['leach_save_images/',xm_str,' nodes alive vs. rounds.png'];
% saveas(third,[figure_name]);
% % 
% % 
% % %====================================================================
% fourth =figure;
% 
% plot(cumsum(PACKETS_TO_BS),100*ones(1,rmax)-DEAD)
% 
% title_label = 'Plot for Network Diameter = ';
% title([title_label xm_str])
% 
% xlabel('Number of data items received at BS (Fig 8(b) TWC)')
% ylabel('Number of nodes alive')
% 
% figure_name = ['leach_save_images/',xm_str,' nodes alive vs. data items at BS.png'];
% saveas(fourth,[figure_name]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   STATISTICS    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                     %
%  DEAD  : a rmax x 1 array of number of dead nodes/round 
%  CLUSTERHS : a rmax x 1 array of number of Cluster Heads/round
%  PACKETS_TO_BS : a rmax x 1 array of number packets send to Base Station/round
%  PACKETS_TO_CH : a rmax x 1 array of number of packets send to ClusterHeads/round
%  first_dead: the round where the first node died                   
%                                                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


minIdx = 5000;
for i=1:length(container_for_plot)
    minIdx = min(container_for_plot(i).lastIdx, minIdx);
end
minIdx

energycontainer = [];
for i=1:length(container_for_plot)
    array = container_for_plot(i).second;
    energycontainer = [energycontainer array(minIdx)];
end

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
%------------------------------------------------------------------------
zero = figure;
for i=1:length(energyplt)
    plot(r_label,energyplt(i).first,line_types(i,:),'LineWidth',1.5);
    hold on
end
legend(legend_f1,legend_f2,legend_f3,legend_f4,legend_f5,legend_f6,legend_f7,legend_f8,legend_f9,legend_f10,legend_f11,'Location','southeast');
title_label = 'Plot for Different Network Size values';
title([title_label])

xlabel('Time(rounds)')
ylabel('Accumulative consumed energy')

figure_name = 'leach_save_images/Accumulative consumed energy vs. rounds.png';
saveas(zero,[figure_name]);

%------------------------f------------------------
maximum = max(energycontainer);
normalised = energycontainer/ maximum;

zero =figure;
plot(xm_set,normalised,'LineWidth',2);
title_label = 'Normalised Energy Dissipation vs. Different Network Size';
title([title_label])

xlabel('Different NetWork Diameter')
ylabel('Normalised Energy Dissipation')

figure_name = 'leach_save_images/f.png';
saveas(zero,[figure_name]);

%---------------------------------------------------

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

figure_name = 'leach_save_images/data items at BS vs. rounds.png';
saveas(first,[figure_name]);

%------------------------------------------------------------------------
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

figure_name = 'leach_save_images/data items at BS vs. energy.png';
saveas(second,[figure_name]);

%------------------------------------------------------------------------
third = figure;
for i=1:length(container_for_plot)
    plot(container_for_plot(i).third,line_types(i,:),'LineWidth',1.5);
    hold on
end
legend(legend_f1,legend_f2,legend_f3,legend_f4,legend_f5,legend_f6,legend_f7,legend_f8,legend_f9,legend_f10,legend_f11,'Location','southwest');
title_label = 'Plot for Different Network Size values';
title([title_label])

xlabel('Time(rounds) (Fig 8(a) TWC)')
ylabel('Number of nodes alive')

figure_name = 'leach_save_images/nodes alive vs. rounds.png';
saveas(third,[figure_name]);

%------------------------------------------------------------------------
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

figure_name = 'leach_save_images/nodes alive vs. data items at BS.png';
saveas(fourth,[figure_name]);


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
legend('10m','30m','50m','70m','90m','110m','130m','150m','170m','190m','200m');
title_label='The round when specific percentage of nodes die for different diameters';
title([title_label]);
xlabel('Amount of deaths');
ylabel('Number of rounds');

figure_name = ['leach_save_images/death round.png'];
saveas(barfig,[figure_name]);