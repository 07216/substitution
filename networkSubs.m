% This function solves the second stage recourse problem. 
% Input: y ordered up to level. d is the realized demand. h is the per unit
% holding cost. p is the per unit shortage cost. hp and pp is the
% equivalent holding and shortage cost when we transform to a problem with
% zero substitution cost.
% s is the substitution cost matrix. N is the problem size. 
% Output: w is the substitution quantity matrix. up is the leftover
% inventory. um is the shortage inventory level. L is the inventory cost.

function [ w,up,um,L ] = networkSubs( y,d,h,p, hp,pp, s,N )

w=zeros(N,N); % matrix of substitution quantity
up=zeros(1,N); % leftover quantity
um=zeros(1,N); % shortage quantity

for i=1:N
    w(i,i)=min([d(i),y(i)]);
    up(i)=max([0, y(i)-w(i,i)]);
    um(i)=max([0,d(i)-w(i,i)]);
end

for j=1:N
    if um(j)>0
        for i=(j-1):-1:1
            if -hp(i)<=pp(j) && up(i)>0
                if up(i)<um(j) % too many shortage
                    w(i,j)=up(i);
                    up(i)=0;
                    um(j)=um(j)-w(i,j);
                end
                if up(i)>=um(j) % enough for shortage
                    w(i,j)=um(j);
                    um(j)=0;
                    up(i)=up(i)-w(i,j);
                end
            end
        end
    end
end

% inventory cost given order up to level y and realized demand d
L=0;
for i=1:N
    L=L+h(i)*up(i)+p(i)*um(i);
    for j=i:N
        L=L+s(i,j)*w(i,j);
    end
end


end

