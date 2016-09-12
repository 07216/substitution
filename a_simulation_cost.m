% given the order up to level y, use simulation to obtain the cost. In the
% 2nd stage we use greedy algorithm to obtain substitution quantities


function [ simCostOpt ] = a_simulation_cost(y_LDR, dist_flag,cvk,cvd,L_capacity, U_capacity, L_demand, U_demand,x,c,s,h,p,hp,pp,N )

CapNum=100; % number of scenarios for capacity realization
DeNum=100;     % number of scenarios for demand realization
simCost=zeros(CapNum,DeNum);

if dist_flag==1
    for kk=1:CapNum
        sim_capacity=L_capacity+(U_capacity-L_capacity)*rand; % generate capacity for each product
        sim_y=min(y_LDR,x+sim_capacity.');
        for j=1:DeNum
            sim_d=L_demand+(U_demand-L_demand)*rand;
            [ sim_w,sim_up,sim_um,inventory_cost ] = networkSubs( sim_y,sim_d,h,p, hp,pp, s,N );
            simCost(kk,j)=inventory_cost+c*(sim_y-x);
        end
    end
elseif dist_flag==2
    for kk=1:CapNum
    for i=1:N
    sim_capacity(i)=trandn_general( (L_capacity(i)+U_capacity(i))/2, cvk(i), L_capacity(i), U_capacity(i), 1 );
    end
    sim_y=min(y_LDR,x+sim_capacity.');
    for j=1:DeNum
        for i=1:N
        sim_d(i)=trandn_general( (L_demand(i)+U_demand(i))/2, cvd(i), L_demand(i), U_demand(i), 1 );
        end
        [ sim_w,sim_up,sim_um,inventory_cost ] = networkSubs( sim_y,sim_d,h,p, hp,pp, s,N );
        simCost(kk,j)=inventory_cost+c*(sim_y-x);
    end
    end
end
simCostOpt=mean(mean(simCost));

end

