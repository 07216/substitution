

clear all 

N=3;

eta=0.1:0.1:0.5;
tau=0.5:0.1:0.9;
cvk=0.1:0.2:0.5;
cvd=0.1:0.2:0.5;

table=allcomb(eta,tau,cvk,cvd);

KK=size(table,1);

simCostOpt=zeros(1,KK);
opt_cost=zeros(1,KK);
T_LDR=zeros(1,KK);
T_opt=zeros(1,KK);
ratio=zeros(1,KK);

parfor ipar=1:KK
    eta=table(ipar,1);
    tau=table(ipar,2);
    cvk=table(ipar,3)*ones(1,N);
    cvd=table(ipar,4)*ones(1,N);
    
% intial inventory level
x=zeros(N,1);

% unit ordering cost
c=zeros(1,N);
% eta=0.1;
for i=1:N
c(i)=1+eta*(N-i);
end

s=zeros(N,N); % substitution cost matrix, sij: use product i to satisfy demand j
T=0.5;
T1=0.5;
T2=0.5;
for i=1:N
    for j=(i+1):N
        s(i,j)=T1*c(i)-T2*c(j);
    end
end

% tau is critial ratio(p-c)/(p+h)
% tau=0.7;

% holding cost, if negative, means salvage value
h=-0.5*c;

p=(h*tau+c)/(1-tau); % shortage cost, calculated based on critial ratio

hp=h-T1*c; % h'_i=h_i-alpha_i; need to be increasing in i
pp=p-T2*c; % p'_j=p_j-beta_j; need to be decreasing in j

for i=1:(N-1)
    if hp(i)>hp(i+1)
        disp('hp is not increasing');
    end
    if pp(i)<pp(i+1)
        disp('pp is not decreasing');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% both capacity and demand are uniform
L_capacity=50*ones(1,N);
U_capacity=150*ones(1,N);

L_demand=50*ones(1,N);
U_demand=150*ones(1,N);

% coefficient of variation 

dist_flag=2; 

trun_k=10;
trun_d=10;

[y_LDR,opt_LDR, T_LDR(ipar)]=a_compute_pieceLDR(dist_flag,cvk,cvd,L_capacity, U_capacity, L_demand, U_demand, trun_k, trun_d,x,c,s,h,p);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Given the optimal order up to level, use simulation to obtain the
% expected cost. In the substitution stage, we use greedy algorithm.
[ simCostOpt(ipar) ] = a_simulation_cost(y_LDR, dist_flag,cvk,cvd,L_capacity, U_capacity, L_demand, U_demand,x,c,s,h,p,hp,pp,N );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute the true optimal cost
Nk=100;
Nd=100;
[ y_opt, opt_cost(ipar), T_opt(ipar) ] = a_compute_true_optimal(Nk,Nd,dist_flag,cvk,cvd,L_capacity, U_capacity, L_demand, U_demand, x,c,s,h,p );

ratio(ipar)=(simCostOpt(ipar)-opt_cost(ipar))/opt_cost(ipar);

end


mean_T_LDR=mean(T_LDR);
mean_T_opt=mean(T_opt);
mean_ratio=mean(ratio);

fileID1 = fopen('result1.txt','w');
fprintf(fileID1, '%f \n\n',N);
fprintf(fileID1,'%f \n',mean_T_LDR);
fprintf(fileID1,'%f \n',mean_T_opt);
fprintf(fileID1,'%f \n',mean_ratio);

fclose(fileID1);


fileID2 = fopen('result2.txt','w');

for i=1:KK
fprintf(fileID1,'%f \n',i);
fprintf(fileID1,'(%f \n',table(i,:));
fprintf(fileID1,'%f ',T_LDR(i));
fprintf(fileID1,'%f ',T_opt(i));
fprintf(fileID1,'%f \n\n',ratio(i));
end
fclose(fileID2);

save('record1','ratio');

% for ipar=1:KK
% fprintf(fileID,'%f \n',ipar);
% fprintf(fileID,'%f \n',ratio(ipar));
% fprintf(fileID,'%f \n',simCostOpt(ipar));
% fprintf(fileID,'%f \n',opt_cost(ipar));
% fprintf(fileID,'%f \n',T_LDR(ipar));
% fprintf(fileID,'%f \n\n',T_opt(ipar));
% 
% end




