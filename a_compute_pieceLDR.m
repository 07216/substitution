
% dist_flag indicate what distribution is used
% 1: uniform, 2: normal

function [y,opt_LDR, T]=a_compute_pieceLDR(dist_flag, cvk,cvd,L_capacity, U_capacity, L_demand, U_demand, trun_r, trun_rd,x,c,s,h,p)

% number of products
N=size(c,2);


mean_d=(L_demand'+U_demand')/2;

% truncation point of capacity
trun_z=zeros(N,trun_r-1); % each row store the truncation points of capacites
% the number of truncation points is turn_r-1

% truncation point of demands

if trun_rd>=2
    trun_zd=zeros(N,trun_rd-1);
else
    trun_zd=0; % do not truncate demand
end

% obtain all truncation points
for i=1:N
    tempGap=(U_capacity(i)-L_capacity(i))/trun_r;
    for j=1:(trun_r-1)
        trun_z(i,j)=L_capacity(i)+tempGap*j;
    end
    
    if trun_rd>=2
        tempGapD=(U_demand(i)-L_demand(i))/trun_rd;
        for j=1:(trun_rd-1)
            trun_zd(i,j)=L_demand(i)+tempGapD*j;
        end
    end
end

% find the polyhedron bound parameters W and h

num_rv2=N*trun_r; % number of lifted variables in the second period
num_rv3=N*trun_rd;
num_rv=1+num_rv2+num_rv3; % number of lifted random variables

[W,H, num_support]=find_support(L_capacity, U_capacity, L_demand, U_demand, trun_r, trun_z, trun_rd, trun_zd);

% Obtain the 2nd order moment matrix via simulation

if dist_flag==1
    [E] = simulation_piecewise_mean_uniform( N, L_capacity, U_capacity, trun_r, trun_z, L_demand, U_demand, trun_rd, trun_zd );
elseif dist_flag==2
    [E] = simulation_piecewise_mean_normal( N, L_capacity, U_capacity, trun_r, trun_z, L_demand, U_demand, trun_rd, trun_zd, cvk, cvd );
end

E2=E(2:num_rv2+1);
E3=E(num_rv2+2:end);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization for LDR
c2=c+h;
C2=[c2.'  zeros(N,N)];  % size N*(N+1)

c3s=zeros(N,N); % store the equivalent cost for w_ij
for i=1:N
    for j=i:N
        c3s(i,j)=s(i,j)-h(i)-p(j);
    end
end

c3=[];
for i=1:N
    c3=[c3, c3s(i,:)];
end


C3=[c3.'  zeros(N^2, 2*N)]; % size (N^2)*(2N+1)

P1=[1 zeros(1,2*N)]; % size 1*(2N+1)
P2=[eye(N+1) zeros(N+1,N)]; % size (N+1)*(2N+1) Notice that P3=I, (2N+1)*(2N+1)

% obtain p'_t
P1p=[1 zeros(1,num_rv-1)];
P2p=[eye(num_rv2+1) zeros(num_rv2+1,num_rv3)];

B=[zeros(N,N+1) eye(N)]; % size N*(2N+1)


% get matrix for the constraints

% to put all decisions in a column vector, for w_ij, we make it [w_11 w_12
% w_13, ... w_n1, ...w_nn]^T.
AAP=zeros(N,N^2);
AAM=zeros(N,N^2);
for i=1:N
    for j=((i-1)*N+i):(i*N)
        AAP(i,j)=1;
    end
end
for i=1:N
    for j=i:N:((i-1)*N+i)
        AAM(i,j)=1;
    end
end

tempA33=[-AAP; -AAM]; % 2N*(N^2)

tempA32=[eye(N)    % 2N*N
     zeros(N,N)];
 
tempB3=[zeros(N,2*N+1)            % 2N*(2N+1)
    zeros(N,N+1) -eye(N)];

A33=[tempA33;eye(N^2)]; % (N^2+2N)*(N^2)
A32=[tempA32;zeros(N^2,N)]; % (N^2+2N)*N
B3=[tempB3;zeros(N^2,2*N+1)];

A22=[eye(N)
    -eye(N)]; % 2N*N

num_constraint2=size(A22,1); % number of constraints for the 2nd stage
num_constraint3=size(A33,1); % number of constraints for the 3rd stage

B2=[x zeros(N,N)
    -x -eye(N)]; % 2N*(N+1)

% Obtain the retraction matrix
R(1,1)=1;
for i=1:N
    temp=ones(1,trun_r);
    R=blkdiag(R,temp);
end

for i=1:N
    temp1=ones(1,trun_rd);
    R=blkdiag(R,temp1);
end

%  store the number of decision variables
dim_dv2=N;
dim_dv3=2*N+N^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
cvx_begin quiet

variable Xw(N^2, num_rv) % this is X3
variable Xv(dim_dv2,num_rv2+1)  % this is X2
variable Lambda2(num_constraint2, num_support)
variable Lambda3(num_constraint3, num_support)

% require that v_i only depend on k_i, and v_i is increasing
Xv(1,2:(trun_r+1))>=0;
Xv(1,(trun_r+2):end)==0;
for i=2:N
    Xv(i,((i-1)*trun_r+2):(i*trun_r+1))>=0;
    Xv(i,2:((i-1)*trun_r+1))==0;
    Xv(i,(i*trun_r+2):end)==0;
end

A22*Xv*P2p-B2*P2*R==Lambda2*W;
A32*Xv*P2p+A33*Xw-B3*R==Lambda3*W;

Lambda2*H>=0;
Lambda2>=0;
Lambda3*H>=0;
Lambda3>=0;

minimize c2*Xv*P2p*E+c3*Xw*E-c*x+p*mean_d

cvx_end

% disp('primal')
% cvx_optval
XXV=full(Xv);
% XXw=full(Xw)
T=toc;

% % obtain the order up to level y
y=XXV(:,1);
for i=1:N
    tempStep=(U_capacity(i)-L_capacity(i))/trun_r;
    if XXV(i,1)>=1
        disp('caution: nonzero intercept')
    end
    if XXV(i,1+(i-1)*trun_r+1)==0
        disp('caution: lower bound is too high')
    else
        y(i)=y(i)+L_capacity(i);
        for j=1:trun_r
            y(i)=y(i)+tempStep*XXV(i,1+(i-1)*trun_r+j);
        end
    end
end

opt_LDR=cvx_optval;

end



