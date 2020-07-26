function zffeval_2(foldername,GCI_Method)
%
% Given a foldername, does evaluation for given epoch
% locations
% addpath(genpath(foldername));

eggs=dir(strcat(foldername,'/*egg*'));
disp('Evaluating the Folder for Epoch Extraction');

IDR=[]; MR=[]; FAR=[]; IDA=[];
for i=1:length(eggs)
    try
        eggs(i).name
        
        [s, fs]=wavread(fullfile(foldername,strrep(eggs(i).name,'-egg','')));
        
        s=resample(s,8000,fs);fs=8000;
        
        
        [polarity] = RESKEW_PolarityDetection(s,fs);
        
        s=s*polarity;
        
        s=s-mean(s); s=s./max(abs(s)*1.01);
        
        [x, fs]=wavread(fullfile(foldername,eggs(i).name));
        
        x=resample(x,8000,fs);fs=8000;egg=x*polarity;
        
        speech1=s(:,1);  % Speech Signal

        time=(1:length(speech1))*1000/fs;    % in milliseconds
        
        degg=diff(egg); degg(end+1)=degg(end);  %dEGG
        
        degg1=degg./max(abs(degg)); 
        
        pind=[];
        
        vad=speech1.^2; vad=smooth(vad,800); vad=vad./max(vad); vad(vad<mean(vad)/20)=0; vad(vad>=mean(vad)/20)=1; 
        vad(1)=0; vad(end)=0; vad=durational_constraint(vad,100*fs/1000);
        
%         vad=getVad_SFF(speech1,fs);
    
        pind=find((diff(vad)==1));
        
        nind=find((diff(vad)==-1));
        
        if(isempty(pind))
            pind=1;
            nind=length(speech1);
        end
        
        for iii=1:length(pind)
            
            speech=speech1(pind(iii):nind(iii));
            degg=degg1(pind(iii):nind(iii));
            
            [~,locs] = findpeaks(-degg,'MINPEAKHEIGHT',0.15,'MINPEAKDISTANCE',0.0025*fs);
            epochs_gt=zeros(size(degg));
            epochs_gt(locs)=1;  % Epochs Ground Truth
            
            if(GCI_Method==1) %% ZFF Method
                [WinLen]=xcorrWinLen(speech,fs);    % Window length for trend removal
                [zfSig]=zeroFreqFilter(speech,fs,WinLen); % ZFF Signal
                [epochs]=epochExtract(zfSig);  % Epoch Extraction
            else if(GCI_Method==2) %% Dypsa Method
                    [gci,goi]=dypsa(s,fs);
                    epochs=s-s;epochs(gci)=1;
                else if(GCI_Method==3) %% SEDREAMS_GCIDetection
                        [WinLen]=xcorrWinLen(speech,fs);
                        [gci,goi] = SEDREAMS_GCIDetection(s,fs,round(1000/WinLen));
                        epochs=s-s;epochs(gci)=1;
                    else if(GCI_Method==4) %% yaga
                            [gci,goi]=yaga(s,fs);
                            epochs=s-s;epochs(gci)=1;
                        else
                            [WinLen]=xcorrWinLen(speech,fs);
                            gci=se_vq(speech,fs,round(1000/WinLen));
                            epochs=s-s;epochs(round(gci*fs))=1;
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
            dev=[]; hit=0; miss=0; fa=0;
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
                elseif sum(count)==0
                    miss=miss+1;
                elseif sum(count)>1
                    fa = fa+1;  % False alarm cycles
                end
                
            end
            
            IDR(end+1) = hit*100/length(larcy)
            MR(end+1) = miss*100/length(larcy);
            FAR(end+1) = fa*100/length(larcy);
            IDA(end+1) = std(dev)*1000/fs;
            
            IDR(IDR==0)=[];
            IDR(isnan(IDR)) = [];
            MR(isnan(MR)) = [];
            FAR(isnan(FAR)) = [];
            IDA(isnan(IDA)) = [];
            
            % % %     if(hit*100/length(larcy)<70)
            % % %         delete(fullfile(foldername,strrep(eggs(i).name,'-egg','')))
            % % %         delete(fullfile(foldername,eggs(i).name))
            % % %     end
            
        end
    catch
        
        continue
    end
    
    
end

IDR(isnan(IDR)) = [];
MR(isnan(MR)) = [];
FAR(isnan(FAR)) = [];
IDA(isnan(IDA)) = [];

disp(['Identification Rate:' num2str(mean(IDR)) '%']);
disp(['Miss Rate:' num2str(mean(MR)) '%']);
disp(['False Alarm Rate:' num2str(mean(FAR)) '%']);
disp(['Identification Accuracy:' num2str(mean(IDA)) 'ms']);

end

