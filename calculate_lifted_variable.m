% given a realization r of a random variable, the number of truncation pieces
% trun_r, and truncation points trun_z, calculate the realization of lifted
% variables


function xi=calculate_lifted_variable(r,trun_r,trun_z)
xi=zeros(1,trun_r);
xi(1)=min(r,trun_z(1));
xi(trun_r)=max(r-trun_z(trun_r-1),0);
for jj=2:(trun_r-1)
    xi(jj)=max(min(r,trun_z(jj))-trun_z(jj-1),0);
end


end

