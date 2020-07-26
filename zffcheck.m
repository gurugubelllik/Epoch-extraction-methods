function [zfSig,gci,degg] = zffcheck(filename)

% Input is filename
% Output is plots asked in A2 Q2

[s,fs]=wavread(filename);   % Reading wave file
speech=s(:,1);  % Speech Signal
egg=s(:,2); % EGG signal
time=(1:length(speech))*1000/fs;    % in milliseconds

[WinLen]=xcorrWinLen(speech,fs);    % Window length for trend removal
[zfSig]=zeroFreqFilter(speech,fs,WinLen); % ZFF Signal
[gci]=epochExtract(zfSig);  % Epoch Extraction
degg = diff(egg); degg(end+1)=degg(end);    %dEGG signal

% Plotting 
h=figure;
ax(1) = subplot(411); plot(time, speech./max(speech)); xlabel('time(ms)'); ylabel('Amplitude'); title('Speech Signal');
ax(2) = subplot(412); plot(time, zfSig./max(zfSig)); xlabel('time(ms)'); ylabel('Amplitude'); title('ZFF Signal');
ax(3) = subplot(413); stem(time, gci,'Marker','none'); xlabel('time(ms)'); ylabel('Amplitude'); title('Epochs from ZFF Signal'); ylim([0 1.5]);
ax(4) = subplot(414); plot(time, degg./max(degg)); xlabel('time(ms)'); ylabel('Amplitude'); title('Differenced EGG Signal');
linkaxes(ax,'x');
set(findall(h, '-property', 'FontSize'), 'FontSize', 12, 'fontWeight', 'bold');

end

