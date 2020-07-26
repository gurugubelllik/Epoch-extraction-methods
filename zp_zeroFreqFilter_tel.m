function [zp_zfSig]=zp_zeroFreqFilter_tel(wav,fs,r,winLen,filter_order)
% clear zp_zfSig;

dwav=(wav);
dwav=diff(wav); dwav(end+1)=dwav(end);
dwav=dwav/max(abs(dwav));

N=length(dwav);

if(filter_order==2)
    
    dwav=(wav);
    dwav=dwav/max(abs(dwav));
    
    if(r==1)
        zp_zfSig=filter([0 -1],conv([1 -r],[r -1]),dwav);
    else
        zp_zfSig=filtfilt(1,[1 -r],dwav);
    end
    
    winLength=winLen;
    zp_zfSig=remTrend(zp_zfSig,winLength*fs/1000);
    zp_zfSig=remTrend(zp_zfSig,winLength*fs/1000);
    
else if(filter_order==4)
        
        if(r==1)
            zp_zfSig=filter([0 -1],conv([1 -r],[r -1]),dwav);
            zp_zfSig=filter([0 -1],conv([1 -r],[r -1]),zp_zfSig);
        else
            zp_zfSig=filtfilt(1,[1 -r],dwav);
            zp_zfSig=filtfilt(1,[1 -r],zp_zfSig);
        end
        
        winLength=winLen;
        
        zp_zfSig=remTrend(zp_zfSig,winLength*fs/1000);
        zp_zfSig=remTrend(zp_zfSig,winLength*fs/1000);
        zp_zfSig=remTrend(zp_zfSig,winLength*fs/1000);
        zp_zfSig=remTrend(zp_zfSig,winLength*fs/1000);
    else if(filter_order == 6)
            
            if(r==1)
                zp_zfSig=filter([0 -1],conv([1 -r],[r -1]),dwav);
                zp_zfSig=filter([0 -1],conv([1 -r],[r -1]),zp_zfSig);
                zp_zfSig=filter([0 -1],conv([1 -r],[r -1]),zp_zfSig);
            else
                zp_zfSig=filtfilt(1,[1 -r],dwav);
                zp_zfSig=filtfilt(1,[1 -r],zp_zfSig);
                zp_zfSig=filtfilt(1,[1 -r],zp_zfSig);
                
            end
            
            winLength=winLen;
            
            zp_zfSig=remTrend(zp_zfSig,winLength*fs/1000);
            zp_zfSig=remTrend(zp_zfSig,winLength*fs/1000);
            zp_zfSig=remTrend(zp_zfSig,winLength*fs/1000);
            zp_zfSig=remTrend(zp_zfSig,winLength*fs/1000);
            zp_zfSig=remTrend(zp_zfSig,winLength*fs/1000);
            zp_zfSig=remTrend(zp_zfSig,winLength*fs/1000);
        else if(filter_order ==8)
                
                if(r==1)
                    zp_zfSig=filter([0 -1],conv([1 -r],[r -1]),dwav);
                    zp_zfSig=filter([0 -1],conv([1 -r],[r -1]),zp_zfSig);
                    zp_zfSig=filter([0 -1],conv([1 -r],[r -1]),zp_zfSig);
                    zp_zfSig=filter([0 -1],conv([1 -r],[r -1]),zp_zfSig);
                else
                    zp_zfSig=filtfilt(1,[1 -r],dwav);
                    zp_zfSig=filtfilt(1,[1 -r],zp_zfSig);
                    zp_zfSig=filtfilt(1,[1 -r],zp_zfSig);
                    zp_zfSig=filtfilt(1,[1 -r],zp_zfSig);
                end
                
                winLength=winLen;
                
                zp_zfSig=remTrend(zp_zfSig,winLength*fs/1000);
                zp_zfSig=remTrend(zp_zfSig,winLength*fs/1000);
                zp_zfSig=remTrend(zp_zfSig,winLength*fs/1000);
                zp_zfSig=remTrend(zp_zfSig,winLength*fs/1000);
                zp_zfSig=remTrend(zp_zfSig,winLength*fs/1000);
                zp_zfSig=remTrend(zp_zfSig,winLength*fs/1000);
                zp_zfSig=remTrend(zp_zfSig,winLength*fs/1000);
                zp_zfSig=remTrend(zp_zfSig,winLength*fs/1000);
            else if(filter_order == 10)
                    if(r==1)
                        zp_zfSig=filter([0 -1],conv([1 -r],[r -1]),dwav);
                        zp_zfSig=filter([0 -1],conv([1 -r],[r -1]),zp_zfSig);
                        zp_zfSig=filter([0 -1],conv([1 -r],[r -1]),zp_zfSig);
                        zp_zfSig=filter([0 -1],conv([1 -r],[r -1]),zp_zfSig);
                        zp_zfSig=filter([0 -1],conv([1 -r],[r -1]),zp_zfSig);
                    else
                        zp_zfSig=filtfilt(1,[1 -r],dwav);
                        zp_zfSig=filtfilt(1,[1 -r],zp_zfSig);
                        zp_zfSig=filtfilt(1,[1 -r],zp_zfSig);
                        zp_zfSig=filtfilt(1,[1 -r],zp_zfSig);
                        zp_zfSig=filtfilt(1,[1 -r],zp_zfSig);
                    end
                    
                    winLength=winLen;
                    
                    zp_zfSig=remTrend(zp_zfSig,winLength*fs/1000);
                    zp_zfSig=remTrend(zp_zfSig,winLength*fs/1000);
                    zp_zfSig=remTrend(zp_zfSig,winLength*fs/1000);
                    zp_zfSig=remTrend(zp_zfSig,winLength*fs/1000);
                    zp_zfSig=remTrend(zp_zfSig,winLength*fs/1000);
                    zp_zfSig=remTrend(zp_zfSig,winLength*fs/1000);
                    zp_zfSig=remTrend(zp_zfSig,winLength*fs/1000);
                    zp_zfSig=remTrend(zp_zfSig,winLength*fs/1000);
                    zp_zfSig=remTrend(zp_zfSig,winLength*fs/1000);
                    zp_zfSig=remTrend(zp_zfSig,winLength*fs/1000);
                else %% if(filter_order == 12)
                    
                    if(r==1)
                        zp_zfSig=filter([0 -1],conv([1 -r],[r -1]),dwav);
                        zp_zfSig=filter([0 -1],conv([1 -r],[r -1]),zp_zfSig);
                        zp_zfSig=filter([0 -1],conv([1 -r],[r -1]),zp_zfSig);
                        zp_zfSig=filter([0 -1],conv([1 -r],[r -1]),zp_zfSig);
                        zp_zfSig=filter([0 -1],conv([1 -r],[r -1]),zp_zfSig);
                        zp_zfSig=filter([0 -1],conv([1 -r],[r -1]),zp_zfSig);
                    else
                        zp_zfSig=filtfilt(1,[1 -r],dwav);
                        zp_zfSig=filtfilt(1,[1 -r],zp_zfSig);
                        zp_zfSig=filtfilt(1,[1 -r],zp_zfSig);
                        zp_zfSig=filtfilt(1,[1 -r],zp_zfSig);
                        zp_zfSig=filtfilt(1,[1 -r],zp_zfSig);
                        zp_zfSig=filtfilt(1,[1 -r],zp_zfSig);
                    end
                    
                    winLength=winLen;
                    
                    zp_zfSig=remTrend(zp_zfSig,winLength*fs/1000);
                    zp_zfSig=remTrend(zp_zfSig,winLength*fs/1000);
                    zp_zfSig=remTrend(zp_zfSig,winLength*fs/1000);
                    zp_zfSig=remTrend(zp_zfSig,winLength*fs/1000);
                    zp_zfSig=remTrend(zp_zfSig,winLength*fs/1000);
                    zp_zfSig=remTrend(zp_zfSig,winLength*fs/1000);
                    zp_zfSig=remTrend(zp_zfSig,winLength*fs/1000);
                    zp_zfSig=remTrend(zp_zfSig,winLength*fs/1000);
                    zp_zfSig=remTrend(zp_zfSig,winLength*fs/1000);
                    zp_zfSig=remTrend(zp_zfSig,winLength*fs/1000);
                    zp_zfSig=remTrend(zp_zfSig,winLength*fs/1000);
                    zp_zfSig=remTrend(zp_zfSig,winLength*fs/1000);
                end
            end
        end
        
    end
end


zp_zfSig(N-240*3:N)=0;
zp_zfSig(1:240*3)=0;
% zp_zfSig=smooth(zp_zfSig,16);




