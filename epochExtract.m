function [gci]=epochExtract(zfSig,sign_flag)


s=zfSig>=0;
k=s(2:end)-s(1:end-1);
if (sign_flag == 1)
    gci=find(k>0);
else
    gci=find(k<0);
end
p=zeros(length(zfSig),1);

% gci=gci+3;
p(gci)=1;
gci=p;








% % 
% % s=zfSig;
% % 
% % s=s./max(abs(s));
% % 
% % s(s<0)=0;
% % 
% % [soe,gci] = findpeaks(s,'MINPEAKHEIGHT',0.0001,'MINPEAKDISTANCE',0.002*8000);
% % % 
% % % s=zfSig>=0;
% % % k=s(2:end)-s(1:end-1);
% % % if (sign_flag == 1)
% % %     gci=find(k>0);
% % % else
% % %     gci=find(k<0);
% % % end
% % p=zeros(length(zfSig),1);
% % 
% % % gci=gci+3;
% % p(gci)=1;
% % gci=p;