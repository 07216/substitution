% Based on find_support, we add more random variables d_k1k2
% trun_dex are the truncation points of the added d_k1k2

function [W,H,num_support,trun_dex]=b_find_support(L_capacity, U_capacity, L_demand, U_demand, trun_r, trun_z, trun_rd, trun_zd)

N=size(L_capacity,2);

% store the lower and upper bound of d_klk2
L_d_ex=zeros(N,N);
U_d_ex=zeros(N,N);
trun_dex=zeros(N,N,trun_rd-1); % store the truncation points of d_k1k2

for k1=1:(N-1)
    for k2=(k1+1):N
        
        L_d_ex(k1,k2)=sum(L_demand(k1:k2));
        U_d_ex(k1,k2)=sum(U_demand(k1:k2));
        tempGap=(U_d_ex(k1,k2)-L_d_ex(k1,k2))/trun_rd;
        for i=1:(trun_rd-1)
            trun_dex(k1,k2,i)=L_d_ex(k1,k2)+tempGap*i;
        end
    end
end


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



for k1=1:(N-1)
    for k2=(k1+1):N
        Vdex=find_lifted_V( L_d_ex(k1,k2),U_d_ex(k1,k2), trun_rd,trun_dex(k1,k2,:) );
        invVdex=inv(Vdex);
        W=blkdiag(W,invVdex(:,2:end));
        H=[H;-invVdex(:,1)];
        
    end
end


num_support=size(W,1); % number of constraints to describe the support
end

