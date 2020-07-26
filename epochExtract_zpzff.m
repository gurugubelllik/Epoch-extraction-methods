function [gci]=epochExtract_zpzff(zfSig)


s=-zfSig;

s=s./max(abs(s));

s(s<0)=0;

% [soe,gci]=findpeaks(s);
% % 
% % % s=filtfilt(blackman(20),1,s);
% % % % % 
[soe,gci] = findpeaks(s,'MINPEAKDISTANCE',0.002*8000);


% [soe,gci] = findpeaks(s,'MINPEAKHEIGHT',0.005,'MINPEAKDISTANCE',0.002*8000);

p=zeros(length(zfSig),1);
p(gci)=1;
gci=p;