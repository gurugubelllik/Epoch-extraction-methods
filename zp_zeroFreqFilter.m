function [zp_zfSig,zp_zff]=zp_zeroFreqFilter(wav,fs,r)
% clear zp_zfSig;
dwav=(wav)-mean(wav);
% dwav=diff(wav); dwav(end+1)=dwav(end);
dwav=dwav/max(abs(dwav));

N=length(dwav);

if(r==1)
    zp_zfSig=filter([0 -1],conv([1 -r],[r -1]),dwav);
    zp_zfSig=filter([0 -1],conv([1 -r],[r -1]),zp_zfSig);
else
    
    zp_zfSig=filtfilt(1,[1 -r],dwav);
    zp_zfSig=filtfilt(1,[1 -r],zp_zfSig);    
    
%     zp_zfSig=filtfilt(1,[1 -r],dwav);
%     zp_zfSig=fliplr(zp_zfSig);
%     zp_zfSig=filtfilt(1,[1 -r],zp_zfSig);
%     zp_zfSig=fliplr(zp_zfSig);
%     zp_zfSig=filtfilt(1,[1 -r],zp_zfSig);
%     zp_zfSig=fliplr(zp_zfSig);
%     zp_zfSig=filtfilt(1,[1 -r],zp_zfSig);
%     zp_zfSig=fliplr(zp_zfSig);
    
end
zp_zff=zp_zfSig;
% max(zp_zfSig)
zp_zfSig=remTrend(zp_zfSig,15*fs/1000);
zp_zfSig=remTrend(zp_zfSig,15*fs/1000);
zp_zfSig=remTrend(zp_zfSig,15*fs/1000);
zp_zfSig=remTrend(zp_zfSig,15*fs/1000);

% zp_zfSig=smooth(zp_zfSig,16);
% zp_zfSig=smooth(zp_zfSig,8);
zp_zfSig(N-240*3:N)=0;
zp_zfSig(1:240*3)=0;
% zp_zfSig=smooth(zp_zfSig,16);

