function [zfSig]=zeroFreqFilter_tel(wav,fs,winLength,filter_order)
dwav=wav;
dwav=diff(wav); dwav(end+1)=dwav(end); % Difference the speech signal
dwav=dwav/max(abs(dwav));
N=length(dwav);



if (filter_order == 2)
    %
    dwav=wav;
    %dwav=diff(wav); dwav(end+1)=dwav(end); % Difference the speech signal
    dwav=dwav/max(abs(dwav));
    zfSig=((((cumsum(cumsum(dwav)))))); % Pass the differenced speech signal twice through zero-frequency resonator..
    
    
    winLength=round(winLength*fs/1000);
    
    zfSig=remTrend(zfSig,winLength);    %Remove the DC offset introduced by zero-frquency filtering
    zfSig=remTrend(zfSig,winLength);
    
else if(filter_order == 4)
        %
        zfSig=((cumsum(cumsum(cumsum(cumsum(dwav)))))); % Pass the differenced speech signal twice through zero-frequency resonator..
        winLength=round(winLength*fs/1000);
        
        zfSig=remTrend(zfSig,winLength);    %Remove the DC offset introduced by zero-frquency filtering
        zfSig=remTrend(zfSig,winLength);
        zfSig=remTrend(zfSig,winLength);
        zfSig=remTrend(zfSig,winLength);
    else if(filter_order == 6)
            %
            zfSig=cumsum(cumsum(cumsum(cumsum(cumsum(cumsum(dwav)))))); % Pass the differenced speech signal twice through zero-frequency resonator..
            winLength=round(winLength*fs/1000);
            
            zfSig=remTrend(zfSig,winLength);    %Remove the DC offset introduced by zero-frquency filtering
            zfSig=remTrend(zfSig,winLength);
            zfSig=remTrend(zfSig,winLength);
            zfSig=remTrend(zfSig,winLength);
            zfSig=remTrend(zfSig,winLength);
            zfSig=remTrend(zfSig,winLength);
            
        else if(filter_order == 8)
                %
                zfSig=cumsum(cumsum(cumsum(cumsum(cumsum(cumsum(cumsum(cumsum(dwav)))))))); % Pass the differenced speech signal twice through zero-frequency resonator..
                winLength=round(winLength*fs/1000);
                
                zfSig=remTrend(zfSig,winLength);    %Remove the DC offset introduced by zero-frquency filtering
                zfSig=remTrend(zfSig,winLength);
                zfSig=remTrend(zfSig,winLength);
                zfSig=remTrend(zfSig,winLength);
                zfSig=remTrend(zfSig,winLength);
                zfSig=remTrend(zfSig,winLength);
                zfSig=remTrend(zfSig,winLength);
                zfSig=remTrend(zfSig,winLength);
                
            else if(filter_order == 10)
                    %
                    zfSig=cumsum(cumsum(cumsum(cumsum(cumsum(cumsum(cumsum(cumsum(cumsum(cumsum(dwav)))))))))); % Pass the differenced speech signal twice through zero-frequency resonator..
                    
                    
                    winLength=round(winLength*fs/1000);
                    
                    zfSig=remTrend(zfSig,winLength);    %Remove the DC offset introduced by zero-frquency filtering
                    zfSig=remTrend(zfSig,winLength);
                    zfSig=remTrend(zfSig,winLength);
                    zfSig=remTrend(zfSig,winLength);
                    zfSig=remTrend(zfSig,winLength);
                    zfSig=remTrend(zfSig,winLength);
                    zfSig=remTrend(zfSig,winLength);
                    zfSig=remTrend(zfSig,winLength);
                    zfSig=remTrend(zfSig,winLength);
                    zfSig=remTrend(zfSig,winLength);
                else %% filter_order == 12
                    %
                    zfSig=cumsum(cumsum(cumsum(cumsum(cumsum(cumsum(cumsum(cumsum(cumsum(cumsum(cumsum(cumsum(dwav)))))))))))); % Pass the differenced speech signal twice through zero-frequency resonator..
                    
                    
                    winLength=round(winLength*fs/1000);
                    
                    zfSig=remTrend(zfSig,winLength);    %Remove the DC offset introduced by zero-frquency filtering
                    zfSig=remTrend(zfSig,winLength);
                    zfSig=remTrend(zfSig,winLength);
                    zfSig=remTrend(zfSig,winLength);
                    zfSig=remTrend(zfSig,winLength);
                    zfSig=remTrend(zfSig,winLength);
                    zfSig=remTrend(zfSig,winLength);    %Remove the DC offset introduced by zero-frquency filtering
                    zfSig=remTrend(zfSig,winLength);
                    zfSig=remTrend(zfSig,winLength);
                    zfSig=remTrend(zfSig,winLength);
                    zfSig=remTrend(zfSig,winLength);
                    zfSig=remTrend(zfSig,winLength);
                end
            end
        end
    end
end

zfSig(N-winLength*3:N)=0;


