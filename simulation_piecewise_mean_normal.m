% use simulation to find the mean for piecewise linear
% decision rule random variables
% everything follows normal distribution, cvk and cvd are the coefficient
% variation of capacity and demand

function [ E] = simulation_piecewise_mean_normal( N, Lk, Uk, trun_r, trun_z, Ld, Ud, trun_rd, trun_zd, cvk, cvd )

meank=(Lk+Uk)/2;
meand=(Ld+Ud)/2;

K=10000; % K is the number of instances
numRow=1+N*trun_rd+N*trun_r; % total number of random variables
MM=zeros(numRow,K);

for k=1:K
    r=zeros(N,1); % capacity
    d=zeros(N,1); % demand
    xi=zeros(N,trun_r); % lifted random variables.
        dd=zeros(N,trun_rd);
    
    for i=1:N
        r(i) = trandn_general( meank(i), cvk(i), Lk(i), Uk(i), 1 ); % generate an instance of random capacity
        % generate instances of lifted random variables
        
        xi(i,:)=calculate_lifted_variable(r(i),trun_r,trun_z(i,:));

        d(i) = trandn_general( meand(i), cvd(i), Ld(i), Ud(i), 1 ); % generate an instance of random demand

        if trun_rd>=2
            dd(i,:)=calculate_lifted_variable(d(i),trun_rd,trun_zd(i,:));
        end
    end
    
    % put everything in a column vector
    % x=zeros(numRow,1);
    x(1)=1;
    
    for i=1:N
        for jj=1:trun_r
            x=[x;xi(i,jj)];
        end
    end
    
    for i=1:N
        if trun_rd<=1
            x=[x;d(i)];
        else
            for jj=1:trun_rd
                x=[x;dd(i,jj)];
            end
        end
    end
    
    MM(:,k)=x;
    clear x
    
end

E=mean(MM,2);


end


