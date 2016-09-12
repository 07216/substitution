% express the support as Wx>=H. num_support is the totoal number of
% constraints.


function [W,H,num_support]=find_support(L_capacity, U_capacity, L_demand, U_demand, trun_r, trun_z, trun_rd, trun_zd)

N=size(L_capacity,2);


W=[1;-1];
H=[1;-1];

for i=1:N


    V  = find_lifted_V( L_capacity(i), U_capacity(i), trun_r,trun_z(i,:) );
    invV=inv(V);
   
    W=blkdiag(W, invV(:,2:end));
    
    H=[H; -invV(:,1)];

end

if trun_rd>=2
    for i=1:N
        Vd = find_lifted_V( L_demand(i), U_demand(i), trun_rd,trun_zd(i,:) );
        invVd=inv(Vd);
        W=blkdiag(W, invVd(:,2:end));
        
        H=[H; -invVd(:,1)];
    end
else
    for i=1:N
        W=blkdiag(W,[1;-1]);
        H=[H;[L_demand(i);-U_demand(i)]];
    end
end

num_support=size(W,1); % number of constraints to describe the support
end

