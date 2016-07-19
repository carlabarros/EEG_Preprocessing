clear all;
clc;

% load('test_sensory_nofilt_ica','EEGseg','chanlocs');
% load('test_sensory_butter_ica','EEGseg','chanlocs');
% load('test_sensory_firfilt_ica','EEGseg','chanlocs');
% load('test_sensory_firfilt01_ica','EEGseg','chanlocs');

% load('test_sensory_nofilt2_ica','EEGseg');
% load('test_sensory_butter2_ica','EEGseg');
% load('test_sensory_firfilt2_ica');

% load('test_sensory_nofilt_bs');
load('test_sensory_butter_bs');
% load('test_sensory_firfilt_bs');
% 

nsubj=2;
nblk=2;

EEGs=cell(nblk,nsubj);

for b=1:nblk
    for s=1:nsubj
        for epo=1:size(EEGseg{b,s}(45,:,:),3);
            EEGs{b,s}(1,:,epo)= bsxfun(@minus,EEGseg{b,s}(45,:,epo),mean(EEGseg{b,s}(45,:,epo),2));
        end
    end
end

% leg_cwt={'cwt_ctr_voi.fig','cwt_exp_voi.fig';'cwt_ctr_bip.fig','cwt_exp_bip.fig'};
% leg_dwt={'dwt_ctr_voi.fig','dwt_exp_voi.fig';'dwt_ctr_bip.fig','dwt_exp_bip.fig'};
% leg_cmor={'cmor_ctr_voi.fig','cmor_exp_voi.fig';'cmor_ctr_bip.fig','cmor_exp_bip.fig'};
% 
% leg1_cwt={'cwt_ctr_v_energy.fig','cwt_exp_v_energy.fig';'cwt_ctr_b_energy.fig','cwt_exp_b_energy.fig'};
% leg1_dwt={'dwt_ctr_v_energy.fig','dwt_exp_v_energy.fig';'dwt_ctr_b_energy.fig','dwt_exp_b_energy.fig'};
% leg1_cmor={'cmor_ctr_v_energy.fig','cmor_exp_v_energy.fig';'cmor_ctr_b_energy.fig','cmor_exp_b_energy.fig'};


%%
%{
%----------CWT------------------------%

ave_powT=cell(nblk,nsubj);%size(chanlocs,2));
ave_plf=cell(nblk,nsubj);%size(chanlocs,2));

scales=6.5:0.395:101.5; %4:0.25:64
coef=single(zeros(100,length(scales),size(EEGseg{1,2}(1,:,:),2)));

for s=1:nsubj
    for b=1:nblk
        for chan=45; %1:length(chanlocs)
            
            for epo=1:size(EEGseg{b,s},3)
                
                fprintf('\n epo number %d',epo);
                
                cfs= cwt(EEGseg{b,s}(chan,:,epo),scales,'morl','scal');
                
                coef(epo,:,:)=cfs;
                
                if epo==1
                    powT_tmp=zeros(size(cfs,1),1024);
                    plf_tmp=zeros(size(cfs,1),1024);
                end
                
                powT_tmp=powT_tmp+(abs(cfs).^2);
                plf_tmp=plf_tmp+(abs(cfs).^2);
            end
            epos=size(EEGseg{b,s},3);
            
            powT_tmp=powT_tmp+(abs(x).^2);
            plf_tmp=plf_tmp+(t./abs(t));
            
            ave_plf{b,s}=abs(plf_tmp./epos);
            ave_powT{b,s}=abs(powT_tmp./epos);
        end
    end
end

%%
intT=1:500;
baseline=700;
srate=512;
time=round((1E-3*(intT(1)+baseline)*srate)+1:round(1E-3*(intT(length(intT))+baseline)*srate));
pfreq = scal2frq(scales,'morl',1/512);

%%{
for s=1:nsubj
    for b=1:nblk
        
        power=ave_powT{b,s};
        
        spec_cwt_filt=sum(power(:,time),2);
                
        figure;
        subplot('position',[0.1 0.2 0.5 0.7])
        imagesc(intT,pfreq,power(:,time),[min(min(power(:,time))) max(max(power(:,time)))]);%
        colorbar;
        colormap jet;
        set(gca,'YDir','normal');
        
        subplot('position',[0.7 0.2 0.2 0.7])
        plot(spec_cwt_filt, pfreq);
        
        savefig(leg_cwt{b,s});
    end
end


for s=1:nsubj
    for b=1:nblk
        
        power=ave_powT{b,s};
        spec_cwt_filt=sum(power(:,time),2);
        figure;
        plot(pfreq,spec_cwt_filt);
        savefig(leg1_cwt{b,s});
    end
end
%}

%%{

%% ----------DWT------------------------%

ave_powT=cell(nblk,nsubj);
ave_plf=cell(nblk,nsubj);

x=single(zeros(8,size(EEGs{1,2}(1,:,:),2)));


for s=1:nsubj
    for b=1:nblk
        for epo=1:1:size(EEGs{b,s},3)
            
            fprintf('\n epo number %d',epo);
            
            [c,l] = wavedec(EEGs{b,s}(1,:,epo),7,'db44');
            
            a7 = wrcoef('a',c,l,'db44',7);
            d7 = wrcoef('d',c,l,'db44',7);
            d6 = wrcoef('d',c,l,'db44',6);
            d5 = wrcoef('d',c,l,'db44',5);
            d4 = wrcoef('d',c,l,'db44',4);
            d3 = wrcoef('d',c,l,'db44',3);
            d2 = wrcoef('d',c,l,'db44',2);
            d1 = wrcoef('d',c,l,'db44',1);
            
            x(1,:)=a7;
            x(2,:)=d7;
            x(3,:)=d6;
            x(4,:)=d5;
            x(5,:)=d4;
            x(6,:)=d3;
            x(7,:)=d2;
            x(8,:)=d1;
            
            if epo==1
                powT_tmp=zeros(size(x,1),1024);
                plf_tmp=zeros(size(x,1),1024);
            end
            
            powT_tmp=powT_tmp+(abs(x).^2);
            plf_tmp=plf_tmp+(x./abs(x));
        end
        epos=size(EEGseg{b,s},3);
        ave_powT{b,s}=abs(powT_tmp./epos);
        ave_plf{b,s}=abs(plf_tmp./epos);
    end
end


%%
intT=-700:1300;
baseline=700;
srate=512;
time=round((1E-3*(intT(1)+baseline)*srate)+1:round(1E-3*(intT(length(intT))+baseline)*srate));

%%{
for s=1:nsubj
    for b=2%1:nblk
        
        power=ave_powT{b,s};
%         power=ave_plf{b,s};
        
        spec_cwt_filt=sum(power(2:8,time),2);
        
        figure;
        subplot('position',[0.1 0.2 0.5 0.7])
        imagesc(intT,2:8,power(2:8,time),[min(min(power(2:8,time))) max(max(power(2:8,time)))]); %
        colorbar;
        colormap jet;
        set(gca,'YDir','normal');
        
        subplot('position',[0.7 0.2 0.2 0.7])
        plot(spec_cwt_filt, 2:8);
        
%         savefig(leg_dwt{b,s});
        
    end
end
%}

% for s=1:nsubj
%     for b=1:nblk
%         
%         power=ave_powT{b,s};
%         %         power=ave_plf{b,s};
%         spec_cwt_filt=sum(power(2:8,time),2);
%         figure;
%         plot(2:8,spec_cwt_filt);
%         savefig(leg1_dwt{b,s});
%         
%     end
% end

%}

%%{

%% MORLET COMPLEX

%Morlet Wavelet
morlf=64:-0.25:4; %4-60Hz step 0.25Hz
morl=cell(1,length(morlf)); %to hold the morlet wavelets for each frequency

% Replicating the parameters used by Roach and Mathalon, 2008
m=4;
c=7;
% nc=(m*c)/(2*pi); %number of cycles
nc=4.46; % m=4 and c=7
srate=512;
step=1/srate;


%Complex Morlet Wavelet definition
for i=1:length(morlf)
    morl{i}=cmorwavf(-nc/morlf(i),nc/morlf(i),(2*nc/morlf(i))*srate,2*(nc/(m*morlf(i)))^2,morlf(i));
end

ave_powT=cell(nblk,nsubj);
ave_plf=cell(nblk,nsubj);

for s=1:nsubj
    for b=1:nblk
        for epo=1:size(EEGs{b,s},3)
            
            for f=1:length(morlf)
                t(f,:)=step*conv(EEGs{b,s}(1,:,epo),morl{f},'same');% gives the absalized complex time-varying energy of each single trial
            end
            
            if epo==1
                powT_tmp=zeros(size(t,1),1024);
                plf_tmp=zeros(size(t,1),1024);
            end
            
            powT_tmp=powT_tmp+(abs(t).^2);
            plf_tmp=plf_tmp+(t./abs(t));
        end
        epos=size(EEGs{b,s},3);
        ave_powT{b,s}=abs(powT_tmp./epos);
        ave_plf{b,s}=abs(plf_tmp./epos);
    end
end

%%
intT=-700:1300;
baseline=700;
srate=512;
time=round((1E-3*(intT(1)+baseline)*srate)+1:round(1E-3*(intT(length(intT))+baseline)*srate));


%%{
for s=1:nsubj
    for b=2%1:nblk
        
        power=ave_powT{b,s};
%         power=ave_plf{b,s};
        
        spec_cwt_filt=sum(power(:,time),2);
                
        figure;
        subplot('position',[0.1 0.2 0.5 0.7])
        imagesc(intT,morlf,power(:,time),[min(min(power(:,time)))  max(max(power(:,time)))]); %
        colorbar;
        colormap jet;
        set(gca,'YDir','normal');
        
        subplot('position',[0.7 0.2 0.2 0.7])
        plot(spec_cwt_filt, morlf);
%         savefig(leg_cmor{b,s});
        
    end
end
%}

% for s=1:nsubj
%     for b=1:nblk
%         
%         power=ave_powT{b,s};
% %         power=ave_plf{b,s};
%         spec_cwt_filt=sum(power(:,time),2);
%         figure;
%         plot(morlf,spec_cwt_filt);
%          savefig(leg1_cmor{b,s});
%     end
% end


%}