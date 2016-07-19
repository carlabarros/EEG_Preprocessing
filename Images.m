%% Window Power mean

baseline=700;
Fs=512;
% intF=209:225;% 4-8Hz
% intF=369:385;
intF=3;
intT=[100 300];

int=round((intT(1)+baseline)*Fs/1000+1:(intT(2)+baseline)*Fs/1000);

% win1=mean(mean(ave_powT{2,1}(intF,int)));
% win2=mean(mean(ave_powT{2,2}(intF,int)));

win1=mean(mean(ave_plf{2,1}(intF,int)));
win2=mean(mean(ave_plf{2,2}(intF,int)));



erpimage(EEGseg{2,1}(45,:,:),[],[-700 1024 512],'Trials',1,1,'erp',2,'cbar');
erpimage(EEGseg{2,2}(45,:,:),[],[-700 1024 512],'Trials',1,1,'erp',2,'cbar');

erpimage(EEGs{2,2}(1,:,:),[],[-700 1024 512],'Trials',1,1,'erp',2,'cbar');

figure;
erpimage(EEGseg{2,1}(1,:,:),[],[-700 1024 512],'Trials',1,1,'erp','cbar');
figure;
erpimage(EEGseg{1,2}(1,:,:),[],[-700 1024 512],'Trials',1,1,'erp','cbar');



rej=cell(2,2);

for i=1:2
    for j=1:2
        
        for m=1:length(trej{i,j})
            
            if trej{i,i}(1,m)==1
                rej{i,j}(1,length(rej{i,j})+1)=m;
            end
        end
    end
end
