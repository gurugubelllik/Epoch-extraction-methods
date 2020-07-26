clc
close all
clear all

addpath(genpath('/home/kalam/Workspace/0_Speech_Analysis_Codes/covarep-master'))
addpath(genpath('/home/kalam/Workspace/0_Speech_Analysis_Codes/GLOAT'))

Path1='/home/kalam/Workspace/wavfiles/GCI_CMU__DATA_2/';
% Path1='/home/kalam/Workspace/7_Prosody_mod/codes/wavs/';


% [IDR1, MR1, FAR1, IDA1, ERR1, ~, Time1]=zffeval_new(Path1,1);
% % % 
% % % [IDR2, MR2, FAR2, IDA2, ERR2, ~, Time2]=zffeval_new(Path1,2);
% % % 
% % % [IDR3, MR3, FAR3, IDA3, ERR3, ~, Time3]=zffeval_new(Path1,3);
% % % 
% % % [IDR5, MR5, FAR5, IDA5, ERR5, ~, Time5]=zffeval_new(Path1,5);
% % % 
% % % [IDR4, MR4, FAR4, IDA4, ERR4, ~, Time4]=zffeval_new(Path1,4);
% % % 
% % % [IDR6, MR6, FAR6, IDA6, ERR6, ~, Time6]=zffeval_new(Path1,6);
% % % 
% % % [IDR7, MR7, FAR7, IDA7, ERR7, ~, Time7]=zffeval_new_tel2(Path1,7);
% % % 
% % % [IDR8, MR8, FAR8, IDA8, ERR8, IDR_25_8, Time8]=zffeval_new_tel2(Path1,8);

%%

Path1='/home/kalam/Workspace/wavfiles/GCI_CMU__DATA_2_Tel_Final/';
% 
% Path1='/home/kalam/Workspace/wavfiles/GCI_CMU__DATA_2_Tel4/';
% 
% Path1='/home/kalam/Workspace/wavfiles/GCI_CMU_Telphonic_NEW/';
% 
[IDR1, MR1, FAR1, IDA1, ERR1, ~, Time1]=zffeval_new_tel2(Path1,1);
% 
% [IDR2, MR2, FAR2, IDA2, ERR2, ~, Time2]=zffeval_new_tel2(Path1,2);
% 
% [IDR3, MR3, FAR3, IDA3, ERR3, ~, Time3]=zffeval_new_tel2(Path1,3);
% 
% [IDR5, MR5, FAR5, IDA5, ERR5, IDR_25_5, Time5]=zffeval_new_tel2(Path1,5);
% 
% [IDR4, MR4, FAR4, IDA4, ERR4, ~, Time4]=zffeval_new_tel2(Path1,4);
% 
% [IDR6, MR6, FAR6, IDA6, ERR6, ~, Time6]=zffeval_new_tel2(Path1,6);
% 
% [IDR7, MR7, FAR7, IDA7, ERR7, ~, Time7]=zffeval_new_tel2(Path1,7);
% 
% [IDR8, MR8, FAR8, IDA8, ERR8, IDR_25_8, Time8]=zffeval_new_tel2(Path1,8);
% 
% [IDR9, MR9, FAR9, IDA9, ERR9, IDR_25_9, Time9]=zffeval_new_tel2(Path1,9);

%%

rmpath(genpath('/home/kalam/Workspace/0_Speech_Analysis_Codes/covarep-master'))
rmpath(genpath('/home/kalam/Workspace/0_Speech_Analysis_Codes/GLOAT'))
%%

% % % figure
% % % a=Time1(1); b=Time1(2); Time1(1)=b;Time1(2)=a;
% % % 
% % % a=Time5(1); b=Time5(2); Time5(1)=b;Time5(2)=a; a= 0.005:0.01:0.1; Time5=Time5+a;
% % % 
% % % 
% % % Time= [Time6' Time2' Time3' Time4' Time1' Time5'];
% % % createfigure(Time)
% legend('SE-VQ', 'Dypsa', 'SEDre','Yaga','ZFF','ZP-ZFF')

clc


%  
disp(['ZFF'])
disp(['  '])
disp(['Identification Rate:' num2str(mean(IDR1)) '%']);
disp(['Miss Rate:' num2str(mean(MR1)) '%']);
disp(['False Alarm Rate:' num2str(mean(FAR1)) '%']);
disp(['Identification Accuracy:' num2str(mean(IDA1)) 'ms']);
disp(['Identification Error:' num2str(mean(ERR1)) 'ms']);
disp(['  '])
disp(['  '])
% % 
% % 
% % disp(['Dypsa'])
% % disp(['  '])
% % disp(['Identification Rate:' num2str(mean(IDR2)) '%']);
% % disp(['Miss Rate:' num2str(mean(MR2)) '%']);
% % disp(['False Alarm Rate:' num2str(mean(FAR2)) '%']);
% % disp(['Identification Accuracy:' num2str(mean(IDA2)) 'ms']);
% % disp(['Identification Error:' num2str(mean(ERR2)) 'ms']);
% % disp(['  '])
% % disp(['  '])
% % 
% % 
% % 
% % disp(['SEDREAMS'])
% % disp(['  '])
% % disp(['Identification Rate:' num2str(mean(IDR3)) '%']);
% % disp(['Miss Rate:' num2str(mean(MR3)) '%']);
% % disp(['False Alarm Rate:' num2str(mean(FAR3)) '%']);
% % disp(['Identification Accuracy:' num2str(mean(IDA3)) 'ms']);
% % disp(['Identification Error:' num2str(mean(ERR3)) 'ms']);
% % disp(['  '])
% % disp(['  '])
% % 
% % 
% % disp(['YAGA'])
% % disp(['  '])
% % disp(['Identification Rate:' num2str(mean(IDR4)) '%']);
% % disp(['Miss Rate:' num2str(mean(MR4)) '%']);
% % disp(['False Alarm Rate:' num2str(mean(FAR4)) '%']);
% % disp(['Identification Accuracy:' num2str(mean(IDA4)) 'ms']);
% % disp(['Identification Error:' num2str(mean(ERR4)) 'ms']);
% % disp(['  '])
% % disp(['  '])
% % 
% % 
% % 
% % disp(['SE_VQ'])
% % disp(['  '])
% % disp(['Identification Rate:' num2str(mean(IDR6)) '%']);
% % disp(['Miss Rate:' num2str(mean(MR6)) '%']);
% % disp(['False Alarm Rate:' num2str(mean(FAR6)) '%']);
% % disp(['Identification Accuracy:' num2str(mean(IDA6)) 'ms']);
% % disp(['Identification Error:' num2str(mean(ERR6)) 'ms']);
% % disp(['  '])
% % disp(['  '])
% % 
% % 
% % disp(['GEFBA'])
% % disp(['  '])
% % disp(['Identification Rate:' num2str(mean(IDR7)) '%']);
% % disp(['Miss Rate:' num2str(mean(MR7)) '%']);
% % disp(['False Alarm Rate:' num2str(mean(FAR7)) '%']);
% % disp(['Identification Accuracy:' num2str(mean(IDA7)) 'ms']);
% % disp(['Identification Error:' num2str(mean(ERR7)) 'ms']);
% % disp(['  '])
% % disp(['  '])


disp(['ZP_ZFF'])
disp(['  '])
disp(['Identification Rate:' num2str(mean(IDR5)) '%']);
disp(['Miss Rate:' num2str(mean(MR5)) '%']);
disp(['False Alarm Rate:' num2str(mean(FAR5)) '%']);
disp(['Identification Accuracy:' num2str(mean(IDA5)) 'ms']);
disp(['Identification Error:' num2str(mean(ERR5)) 'ms']);
disp(['  '])
disp(['  '])
% % 
% % 
% % disp(['CWT_GCI'])
% % disp(['  '])
% % disp(['Identification Rate:' num2str(mean(IDR8)) '%']);
% % disp(['Miss Rate:' num2str(mean(MR8)) '%']);
% % disp(['False Alarm Rate:' num2str(mean(FAR8)) '%']);
% % disp(['Identification Accuracy:' num2str(mean(IDA8)) 'ms']);
% % disp(['Identification Error:' num2str(mean(ERR8)) 'ms']);
% % disp(['  '])
% % disp(['  '])
% % 
% % 
% % 
% % disp(['SPF'])
% % disp(['  '])
% % disp(['Identification Rate:' num2str(mean(IDR9)) '%']);
% % disp(['Miss Rate:' num2str(mean(MR9)) '%']);
% % disp(['False Alarm Rate:' num2str(mean(FAR9)) '%']);
% % disp(['Identification Accuracy:' num2str(mean(IDA9)) 'ms']);
% % disp(['Identification Error:' num2str(mean(ERR9)) 'ms']);
% % disp(['  '])
% % disp(['  '])
% % 


