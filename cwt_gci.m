function evedence=cwt_gci(s)
c=cwt(hilbert(s),2:21,'sym4');

ev=mean(abs(c));

gwin=diff(gausswin(64,2.8571));

evedence=filtfilt(gwin,1,ev);

end