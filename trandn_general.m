% generate N truncated normal random variables with mean=m, coefficient of
% variation=cv=\sigma/\mu, lower bound=l, upper bound=u. Use outside
% function trandn

function [ Z ] = trandn_general( m, cv, l, u, N )

Z=zeros(N,1);

for i=1:N
s=m*cv;
Z(i)=m+s*trandn((l-m)/s,(u-m)/s);
end

end

