function [IDR, MR, FAR, IDA, ERR, IDR_25, Time]=zffeval_new_tel(foldername,GCI_Method)
%
% Given a foldername, does evaluation for given epoch
% locations
% addpath(genpath(foldername));

FileList = dir([foldername '*.wav']);
disp('Evaluating the Folder for Epoch Extraction');


IDR=[]; MR=[]; FAR=[]; IDA=[]; ERR=[];IDR_25=[]; Time=[];
for i=1:length(FileList)
    try
        fname = FileList(i).name;   % Filename
        disp(['File No. ' num2str(i)])
        %     disp(['Filename: ' fname])
        
        [s,fs]=wavread([foldername fname]);   % Reading wave file
        %     length(s)
        s=resample(s,8000,fs);fs=8000;
        %     length(s)
        speech=s(:,1);  % Speech Signal
        
        %         speech=filtfilt(gausswin(40),1,speech);
        
        egg=s(:,2); % EGG signal
        time=(1:length(speech))*1000/fs;    % in milliseconds
        
        degg=diff(egg); degg(end+1)=degg(end);  %dEGG
        degg=degg./max(abs(degg));
        
        
        [~,locs] = findpeaks(-degg,'MINPEAKHEIGHT',0.15,'MINPEAKDISTANCE',0.001*fs);
        
        %         [xc,lag]=xcorr(speech,degg);
        %         [mx,mind]=max(abs(xc));
        %         delay =lag(mind);
        %         locs=locs-delay;
        %
        epochs_gt=zeros(size(degg));
        epochs_gt(locs)=1;  % Epochs Ground Truth
        
        
        
        
        if(GCI_Method==1) %% ZFF Method
            disp(['ZFF'])
            tic
            [WinLen]=xcorrWinLen(speech,fs);    % Window length for trend removal
            [zfSig]=zeroFreqFilter_tel(speech,fs,WinLen); % ZFF Signal
            [epochs]=epochExtract(zfSig,1);  % Epoch Extraction
            time_elapsed=toc
            
        else if(GCI_Method==2) %% Dypsa Method
                disp(['Dypsa'])
                tic
                [gci,goi]=dypsa(s,fs);
                epochs=s-s;epochs(gci)=1;
                time_elapsed=toc
            else if(GCI_Method==3) %% SEDREAMS_GCIDetection
                    disp(['SEDREAMS'])
                    tic
                    wave=speech;
                    F0min=80;
                    F0max=240;
                    % Pitch tracking using a method based on the Summation of the Residual Harmonics
                    [f0,VUVDecisions,SRHVal] = SRH_PitchTracking(wave,fs,F0min,F0max);
                    
                    VUVDecisions2=zeros(1,length(wave));
                    HopSize=round(10/1000*fs);
                    for k=1:length(VUVDecisions)
                        VUVDecisions2((k-1)*HopSize+1:k*HopSize)=VUVDecisions(k);
                    end
                    % Estimation of the mean pitch value
                    f0_tmp=f0.*VUVDecisions;
                    pos= f0_tmp~=0;
                    f0_tmp=f0_tmp(pos);
                    F0mean=mean(f0_tmp);
                    % Oscillating Moment-based Polarity Detection
                    [Polarity] = RESKEW_PolarityDetection(wave,fs);
                    
                    % Speech Event Detection using the Residual Excitation And a Mean-based Signal
                    [gci,MeanBasedSignal] = SEDREAMS_GCIDetection(Polarity*wave,fs,F0mean);
                    epochs=s-s;epochs(gci)=1;
                    time_elapsed=toc
                else if(GCI_Method==4) %% yaga
                        disp(['YAGA'])
                        tic
                        [gci,goi]=yaga(s,fs);
                        epochs=s-s;epochs(gci)=1;
                        time_elapsed=toc
                    else if(GCI_Method==5) %% zp_zff
                            disp(['ZP_ZFF'])
                            %                             tic
                            %                             [WinLen]=xcorrWinLen(speech,fs);    % Window length for trend removal
                            [zfSig]=zp_zeroFreqFilter_tel(speech,fs,0.98,15); % ZFWinLenF Signal
                            [epochs]=epochExtract_zpzff(zfSig);  % Epoch Extraction
                            %                             time_elapsed=toc;
                        else if(GCI_Method==6) %%GEFBA
                                disp(['SE_VQ'])
                                tic
                                [WinLen]=xcorrWinLen(speech,fs);
                                gci=se_vq(speech,fs,round(1000/WinLen));
                                epochs=s-s;epochs(round(gci*fs))=1;
                                time_elapsed=toc
                                
                            else%% se_vq
                                disp(['GEFBA'])
                                [GOIs, gci, ~, ~, ~, ~, ~]=GEFBA(s,fs);
                                epochs=s-s;epochs(gci)=1;
                            end
                        end
                    end
                end
            end
        end
        % Larynx Cycles
        temp=find(epochs_gt==1);
        N=length(temp);
        init = [temp(1:N-2)+temp(2:N-1)]/2;
        fin = [temp(2:N-1)+temp(3:N)]/2;
        larcy=[round(init) round(fin)];   % Larynx Cycles
        
        % Duration Constraint
        temp = fin-init<=(15*fs/1000);  % 15ms
        ind = find(temp==1);
        larcy = larcy(ind,:);
        
        epochlocs=find(epochs==1);
        epochlocs_gt=find(epochs_gt==1);
        dev=[]; hit=0; miss=0; fa=0; hit1=0;
        for j=1:length(larcy)
            
            % Find detected epochs inside larynx cycle i
            count1 = (larcy(j,1)<=epochlocs);
            count2 = (epochlocs<=larcy(j,2));
            count = and(count1, count2);    % Binary vector indicating detected epochs in given larynx cycle
            
            if sum(count)==1
                hit=hit+1;
                
                % Find ground truth epoch near detected epoch
                ind1 = (larcy(j,1)<=epochlocs_gt);
                ind2 = (epochlocs_gt<=larcy(j,2));
                ind = and(ind1, ind2);      % Binary vector indicating epoch in ground truth in given larynx cycle
                
                dev(end+1) = abs(epochlocs_gt(find(ind==1))-epochlocs(find(count==1)));   % Accumulate the deviations
                
                if (abs(dev)<5)
                    hit1=hit1+1;
                end
            elseif sum(count)==0
                miss=miss+1;
            elseif sum(count)>1
                fa = fa+1;  % False alarm cycles
            end
            
        end
        hit*100/length(larcy)
        IDR_25(end+1) = hit1*100/length(larcy);
        IDR(end+1) = hit*100/length(larcy);
        MR(end+1) = miss*100/length(larcy);
        FAR(end+1) = fa*100/length(larcy);
        IDA(end+1) = std(dev)*1000/fs;
        ERR(end+1)= mean(dev);
        Time(end+1)=time_elapsed;
        
        IDR(isnan(IDR)) = [];
        MR(isnan(MR)) = [];
        FAR(isnan(FAR)) = [];
        IDA(isnan(IDA)) = [];
        ERR(isnan(ERR)) = [];
        Time(isnan(IDR)) = [];
        % % %     if(hit*100/length(larcy)<70)
        % % %         delete(fullfile(foldername,strrep(eggs(i).name,'-egg','')))
        % % %         delete(fullfile(foldername,eggs(i).name))
        % % %     end
    catch
        
        continue
    end
end

IDR_25(isnan(IDR_25)) = [];
IDR(isnan(IDR)) = [];
MR(isnan(MR)) = [];
FAR(isnan(FAR)) = [];
IDA(isnan(IDA)) = [];
ERR(isnan(ERR)) = [];
Time(isnan(IDR)) = [];

disp(['Identification Rate:' num2str(mean(IDR)) '%']);
disp(['Miss Rate:' num2str(mean(MR)) '%']);
disp(['False Alarm Rate:' num2str(mean(FAR)) '%']);
disp(['Identification Accuracy:' num2str(mean(IDA)) 'ms']);
disp(['Identification Error:' num2str(mean(ERR)) 'ms']);

end

