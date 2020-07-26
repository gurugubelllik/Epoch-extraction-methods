function wav=pcmread(wavname)

fid=fopen(wavname);
sig=double(int16(fread(fid,inf,'int16')));
wav=sig./max(1.1*abs(sig));
fclose(fid);

end