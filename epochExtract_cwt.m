function [gci]=epochExtract_cwt(evidence)


s=evidence;

s=s./max(abs(s));

% s(s<0)=0;

[~,gci] = findpeaks(s,'MINPEAKDISTANCE',0.002*8000);

p=zeros(length(evidence),1);
p(gci)=1;
gci=p;