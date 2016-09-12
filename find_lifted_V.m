% given the lower bound L, upper bound U, number of truncation pieces
% trun_r, and the truncation points trun_z, find the matrix V in the paper
% Generalized decision rule approximation for
% stochstic Programming via Lifting to construct convex hull of the lifted
% space.


function [ V ] = find_lifted_V( L, U, trun_r,trun_z )
V=zeros(trun_r+1,trun_r+1);
V(1,:)=ones(1,trun_r+1);
V(2,1)=L;
V(2,2:end)=trun_z(1);
V(trun_r+1,trun_r+1)=U-trun_z(trun_r-1);
for i=3:(trun_r)
    for j=i:(trun_r+1)
        V(i,j)=trun_z(i-1)-trun_z(i-2);
    end
end


end

