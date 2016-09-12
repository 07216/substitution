% 
% y1=50;
% y2=50;
% 
% for i=1:100
%     d1(i)=i;
%     d2(i)=2*i;
%     dsum(i)=d1(i)+d2(i);
%     f(i)=min(max(d2(i)-y2,0),max(y1-d1(i),0));
% end
% 
% plot(dsum,f)
% grid on

clear all

fileID = fopen('result1-2.txt','r');

C = textscan(fileID,'%s','Delimiter','\n');

fclose(fileID);

NN=225;
info=zeros(NN,3);
for i=1:NN
    info(i,:)=str2num(C{1}{7*(i-1)+6});
end
ratio_N3=info(:,3);
 max(ratio_N3)
% save('ratio_N3','ratio_N3');
