function [epochs]=SPF_GCI_Extraction(wav,fs)

y=mean(SFF_SPECTRUM(wav,fs,50,20,3600,0.9394)');

y=y(:);

gwin=diff(gausswin(64,2.666));

evidence=filtfilt(gwin,1,y);

s=evidence;

s=s./max(abs(s));

[~,gci] = findpeaks(s,'MINPEAKDISTANCE',0.002*8000);

p=zeros(length(evidence),1);
p(gci)=1;
epochs=p;

end
