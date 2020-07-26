


wav=s; preprocess_flag=1;


if(preprocess_flag==1)
    F_high=fs/2;    
    F0_high=65;    
    No_levels=round(log2(F_high/F0_high));    
    [WPT]=modwpt(wav,'db5',No_levels,'TimeAlign', true);
    
    dwav=sum(WPT(2:10,:));
else
    dwav=(wav);
end

% dwav=diff(wav); dwav(end+1)=dwav(end);

dwav=dwav/max(abs(dwav));