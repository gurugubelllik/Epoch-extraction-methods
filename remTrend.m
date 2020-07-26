	
function [out]=remTrend(sig,winSize)

window=ones(winSize,1);
rm=conv(sig,window);
rm=rm(winSize/2:length(rm)-winSize/2);


norm=conv(ones(length(sig),1),window);
norm=norm(winSize/2:length(norm)-winSize/2);
	
rm=rm./norm;

out=sig-rm;
