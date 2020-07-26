function vnvcheck(filename)
% 
% Given wav file name as input
% generates plots required for A2 Q3

[s,fs] = wavread(filename);
speech=s(:,1); speech = speech./max(abs(speech));
egg=s(:,2);
time = (1:length(speech))*1000/fs;

windows = buffer(speech,20*fs/1000,(20*fs/1000)-1); % 20ms size, 1sample shift

% Calculating Parameters
E=[]; Z=[]; NACC=[]; PEER=[];
for i=1:size(windows,2)

    E(end+1) = sum([windows(:,i)].^2); % Energy
    
    [zcloc,~]=zerocros(windows(:,i),'p');
    Z(end+1) = length(zcloc);    % ZCR

    num = sum(windows(2:end,i).*windows(1:end-1,i));
    denom = sum(windows(:,i).^2);
    NACC(end+1) = num/denom;    % Normalized Auto Corr Coeff
    
    num = sum(windows(2:end,i)-windows(1:end-1,i));
    denom = sum(windows(:,i).^2);
    PEER(end+1) = num/denom; % Pre-emphasized enery ratio
end
E=E./max(abs(E));
Z=Z./max(abs(Z));
NACC=NACC./max(abs(NACC));
% PEER=PEER./max(abs(PEER));

degg=diff(egg); degg(end+1)=degg(end);
degg=degg./max(abs(degg));

% Taking VNV Decisions
E_VNV = (E>0.1);
Z_VNV = (Z<0.2);
NACC_VNV = (NACC>0.3);
PEER_VNV = (PEER<0.1);



% Plotting

h=figure; 
ax(1)=subplot(311); plot(time,speech,'LineWidth',2.0); xlabel('time(ms)');ylabel('Amplitude'); title('Speech Signal');
ax(2)=subplot(312); plot(time,E,'LineWidth',2.0); hold on; plot(time,0.8.*E_VNV,'r','LineWidth',2.0); xlabel('time(ms)');ylabel('Amplitude'); title('Energy');
ax(3)=subplot(313); plot(time,degg,'LineWidth',2.0); xlabel('time(ms)');ylabel('Amplitude'); title('dEGG');
linkaxes(ax,'x');
set(findall(h, '-property', 'FontSize'), 'FontSize', 12, 'fontWeight', 'bold','LineWidth',2.0);
axis tight;

h=figure; 
ax(1)=subplot(311); plot(time,speech,'LineWidth',2.0); xlabel('time(ms)');ylabel('Amplitude'); title('Speech Signal');
ax(2)=subplot(312); plot(time,Z,'LineWidth',2.0); hold on; plot(time,0.8.*Z_VNV,'r','LineWidth',2.0); xlabel('time(ms)');ylabel('Amplitude'); title('ZCR');
ax(3)=subplot(313); plot(time,degg,'LineWidth',2.0); xlabel('time(ms)');ylabel('Amplitude'); title('dEGG');
linkaxes(ax,'x');
set(findall(h, '-property', 'FontSize'), 'FontSize', 12, 'fontWeight', 'bold');
axis tight;

h=figure; 
ax(1)=subplot(311); plot(time,speech,'LineWidth',2.0); xlabel('time(ms)');ylabel('Amplitude'); title('Speech Signal');
ax(2)=subplot(312); plot(time,NACC,'LineWidth',2.0); hold on; plot(time,0.8.*NACC_VNV,'r','LineWidth',2.0); xlabel('time(ms)');ylabel('Amplitude'); title('Normalized Auto-Correlation Coefficient');
ax(3)=subplot(313); plot(time,degg,'LineWidth',2.0); xlabel('time(ms)');ylabel('Amplitude'); title('dEGG');
linkaxes(ax,'x');
set(findall(h, '-property', 'FontSize'), 'FontSize', 12, 'fontWeight', 'bold','LineWidth',2.0);
axis tight;

h=figure; 
ax(1)=subplot(311); plot(time,speech,'LineWidth',2.0); xlabel('time(ms)');ylabel('Amplitude'); title('Speech Signal');
ax(2)=subplot(312); plot(time,PEER,'LineWidth',2.0); ylim([-2 2]); hold on; plot(time,0.8.*PEER_VNV,'r','LineWidth',2.0); xlabel('time(ms)');ylabel('Amplitude'); title('Pre-Emphasized Energy Ratio');
ax(3)=subplot(313); plot(time,degg,'LineWidth',2.0); xlabel('time(ms)');ylabel('Amplitude'); title('dEGG');
linkaxes(ax,'x');
set(findall(h, '-property', 'FontSize'), 'FontSize', 12, 'fontWeight', 'bold','LineWidth',2.0);
% axis tight;