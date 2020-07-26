function [zfSig]=zeroFreqFilter(wav,fs,winLength)
dwav=wav;
dwav=diff(wav); dwav(end+1)=dwav(end); % Difference the speech signal
dwav=dwav/max(abs(dwav));    
N=length(dwav);

	
zfSig=cumsum(cumsum(cumsum(cumsum(dwav)))); % Pass the differenced speech signal twice through zero-frequency resonator..


winLength=round(winLength*fs/1000);

zfSig=remTrend(zfSig,winLength);    %Remove the DC offset introduced by zero-frquency filtering
zfSig=remTrend(zfSig,winLength);
zfSig=remTrend(zfSig,winLength);
zfSig=remTrend(zfSig,winLength);
zfSig(N-winLength*3:N)=0;

	
