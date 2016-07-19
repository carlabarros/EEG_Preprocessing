%%{
clear all;
clc;

eeglab;
load('channels_79.mat', 'channels');
srate=512;

load('cfg');
load('test_sensory_nofilt2_ica','comprej');

nsubj=2;% Control + experimental
nblk=2;% voices+bips
srej=zeros(1,nsubj); % variable for subject rejection
crej1={'GSR1','GSR2','Erg1','Erg2','Resp','Plet','Temp'};
chn_elec={'FP1','AF7','AFz'};
chn_ext={'EXG6','EXG7','EXG8'};
crej={'EXG6','EXG7','EXG8','GSR1','GSR2','Erg1','Erg2','Resp','Plet','Temp'}; % variable for channel rejection
cinterp=cell(1,nsubj); %variable for channel interpolation
autocint=1;

ref=[65 66]; % mastoid indexes for referencing

bounds=[-0.7 1.3]; % interval for averages

allTr={'30','70','60','65','75','20','276','40','296','45','301','55'};

ao_trg_voice={'55'};
ao_trg_bip={'75'};

n_ao=100; %number of AO triggers

segcount=zeros(nsubj,nblk); %variable to save amount segments used for averages
seg_file=cell(nsubj,5); % 5 files

str_suf={'_AMVoice_AM_AO_Voice.bdf','_AMVoice_AM_AO_Bip.bdf'};

%------------------------------------------Channels:
chn=1:size(channels,2);

for ii=1:length(crej)
    for i=1:length(channels)
        if strcmp(channels(i).labels,crej{1,ii})==1
            chn(i)=0;
        end
    end
end
chn=sort(chn);
chn=chn(length(crej)+1:length(chn));

chanlocs=channels(chn(1:length(chn)));

%-----------------------------------------subj 1 and 2:
chn1=1:size(channels,2);
chn2=1:size(chn_elec,2);
chn3=1:size(chn_ext,2);

for ii=1:length(crej1)
    for i=1:length(channels)
        if strcmp(channels(i).labels,crej1{1,ii})==1
            chn1(i)=0;
        end
        for j=1:3
            if strcmp(channels(i).labels,chn_elec{1,j})==1
                chn2(j)=i;
            end
            if strcmp(channels(i).labels,chn_ext{1,j})==1
                chn3(j)=i;
            end
            
        end
    end
end
chn1=sort(chn1);
chn1=chn1(length(crej1)+1:length(chn1));


%---------------------------------------------Variables

EEGseg=cell(nblk,nsubj); %matrix EEG(blocks, subjects)

rmepo2=cell(nsubj);
comprej2=cell(nsubj);
trej2=cell(nsubj,nblk);


ave_aov=cell(1,size(chanlocs,2)-5);% ao_voices
ave_aob=cell(1,size(chanlocs,2)-5);% ao_bips

%}
for strind=1:nsubj %cycle through all the subjects

    if srej(strind)==0 %unless those selected for rejection
        
        fprintf('\n\nNow processing subject %d\n',strind)
        
        if strind==1
            str='C:\Users\Carla\OneDrive\Documentos\PhD Basic Psychology\Courses\Optional I- Técnicas Avançadas de Imagem Médica\Teste_Sensory\Control\C008'; %prepares string with file name
        elseif strind==nsubj
            str='C:\Users\Carla\OneDrive\Documentos\PhD Basic Psychology\Courses\Optional I- Técnicas Avançadas de Imagem Médica\Teste_Sensory\Experimental\C029'; %prepares string with file name
        end
        
        subjdone=0; % variable to infer if the subject has been processed (1), requires interpolation (0) or is to be rejected (-1)
        
        while subjdone==0
            
            subjdone=1;
            
            % variables to hold subject segments for each condition
            eeg1={}; %AO_voices            
            eeg2={}; %AO_bips
            
            strt1=str;
            strt2=str;
            
            if subjdone==1 || length(cinterp{strind})>4 % unless the subjected has been targeted for interpolation or rejection
                
                %----load file-----%
                strt1(length(strt1)+1:length(strt1)+length(str_suf{1}))=str_suf{1};
                
                EEG1 = pop_fileio(strt1);
                
                strt2(length(strt2)+1:length(strt2)+length(str_suf{2}))=str_suf{2};
                
                EEG2 = pop_fileio(strt2);
                
                
                %--- move AO_voices triggers 300ms backward---%
                for i=1:length(EEG1.event)
                    if EEG1.event(i).type==55
                        EEG1.event(i).latency = EEG1.event(i).latency-0.3*srate;
                    end
                end
                
                for i=1:length(EEG2.event)
                    if EEG2.event(i).type==55 || EEG2.event(i).type==311
                        EEG2.event(i).type=75;
                        EEG2.event(i).latency = EEG2.event(i).latency-0.05*srate; %move triggers 50ms backward
                    end
                end
                
                %---------------Merging Datasets----------------------------%
                
                EEG= pop_mergeset(EEG1,EEG2,1);
                                
                %-------------Select Channels --------------------------%
                
                EEG= pop_select(EEG,'channel',chn); %from the 69 channels only the earlier selected are used
                
                EEG.chanlocs=chanlocs'; %save channel location onto variable
                
                %--------------- Re-referencing -------------------------%
                
                EEG= pop_reref(EEG, ref);
                
                %------- Highpass Filter 1.2Hz - FIRFILT ----------%
                
                dev=0.001;
                fs = 512;
                cutoff = 0.1;%1.2Hz
                df=cutoff*1.001;
                
                beta=pop_kaiserbeta(dev);
                m= pop_firwsord('kaiser', fs, df, dev);
                [EEG,comd, b] = pop_firws(EEG,'fcutoff', cutoff, 'ftype', 'highpass', 'wtype', 'kaiser', 'warg', beta, 'forder', m, 'minphase', 0, 'plotfresp',1);

               
                
%                 [EEG,comd, b] = pop_eegfiltzpButter(EEG,0.1,[],4,0); % 4th order high-pass filter - cutoff:0.1Hz
                
                
                % ------- Epoching ---------------------%
                
                EEG = pop_epoch( EEG, allTr, bounds, 'epochinfo', 'yes');

                %%{
                %-------% Automatic artifact/epoch rejection %-------%
                
                [EEG, rmpo]=pop_autorej(EEG,'eegplot','off','electrodes',[1:EEG.nbchan],'nogui','on', 'startprob',10, 'maxrej', 1);
                                
                rmepo2{strind}=rmpo;
                
                %--------% ICA using Extended InfoMax%------------------------%
                EEG = pop_runica(EEG,'icatype','runica','extended',1,'interupt','on');

               
                %-------- Inspect and reject components ------------------%

                EEG=eeg_SASICA(EEG,cfg);
                close all;
             
                comprej2{strind}=EEG.reject.gcompreject;
%                 EEG = pop_subcomp(EEG,find(EEG.reject.gcompreject),0); %1);
                EEG = pop_subcomp(EEG,find(comprej{strind}),0); %1);
                %} 
                
                blk=1;
                
                if isempty(EEG)==0 %prevents errors by testing if the EEG variable isn't empty
                    
                    %........................................................................................................imag_triggers
                    
                    EEG1 = pop_epoch( EEG, ao_trg_voice, bounds, 'newname', 'BDF file epochs', 'epochinfo', 'yes'); %extracts epochs with the standard triggers of this block
                    
                    [EEG1, rejtrsh1] = eegthresh_simple(EEG1,100,[]); %removes high amplitude events
                                        
                    EEG1 = pop_select( EEG1,'nochannel',{'HEOGL' 'HEOGR' 'VEOG'});
                    
                    eeg1{blk}=EEG1.data; %the data is saved
                                        
                    %........................................................................................................am_triggers
                    
                    EEG2 = pop_epoch(EEG, ao_trg_bip , bounds, 'newname', 'BDF file epochs', 'epochinfo', 'yes'); %same procedure but with the deviant epochs
                    
                    [EEG2, rejtrsh2] = eegthresh_simple(EEG2,100,[]);
                    
                    EEG2 = pop_select( EEG2,'nochannel',{'HEOGL' 'HEOGR' 'VEOG'});
                    
                    eeg2{blk}=EEG2.data;
                    
                end
                
                trej{strind,1}=rejtrsh1;
                trej{strind,2}=rejtrsh2;
                
                EEGseg{1,strind}=eeg1{1};
                EEGseg{2,strind}=eeg2{1};
                
            end
        end
        
    end
end

save('test_sensory_firfilt2_ica','-v7.3');

