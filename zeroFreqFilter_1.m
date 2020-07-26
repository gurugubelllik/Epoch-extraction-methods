function [zfSig,zp_zff]=zeroFreqFilter_1(wav,fs,winLength,r)
clear zfSig;
dwav=wav;
% dwav=diff(wav); dwav(end+1)=dwav(end); % Difference the speech signal
dwav=dwav/max(abs(dwav));    
N=length(dwav);

zfSig=filter(1,conv([1 -r],[1 -r]),dwav);
zfSig=filter(1,conv([1 -r],[1 -r]),zfSig);
% zfSig=cumsum(cumsum(cumsum(cumsum(dwav)))); % Pass the differenced speech signal twice through zero-frequency resonator..

max(zfSig)
winLength=round(winLength*fs/1000);

zp_zff=zfSig;

zfSig=remTrend(zfSig,winLength);    %Remove the DC offset introduced by zero-frquency filtering
zfSig=remTrend(zfSig,winLength);
zfSig=remTrend(zfSig,winLength);
zfSig=remTrend(zfSig,winLength);
zfSig(N-winLength*3:N)=0;

	
