
  virtsens = [];
    for k = 1:length(data_appended.trial)
        virtsens.trial{k} = spatialfilter*data_appended.trial{k};
    end


[pxx,f] = pwelch(spatialfilter,[],[],[],256.001);    
figure
semilogy(f, pxx)

[FourRef,Fref]=fftBasic(spatialfilter', round(256.001));
 figure;
 plot(FourRef, Fref)

 
 [pxx,f] = pwelch(virtsens.trial{1,1}',[],[],[],256.001);    
figure
semilogy(f, pxx)

[FourRef,Fref]=fftBasic(virtsens.trial{1,1}, round(256.001));
 figure;
 plot(FourRef, Fref)

 
  
 [pxx,f] = pwelch(virtsens.trial{1,1}',[],[],[],256.001);    
figure
semilogy(f, pxx)

[FourRef,Fref]=fftBasic(virtsens.trial{1,1}, round(256.001));
 figure;
 plot(FourRef, Fref)




ns_avg=ns_avg';

sensordata_all_subj(8)

cfg=[];
avg=ft_timelockanalysis(cfg, sensordata_all_subj(8));

avg=squeeze(mean(avg_data_dislike.trial));
sRate=256;    
[FourRef,Fref]=fftBasic(sensordata_all_subj(8).trial{1,1},round(sRate));
figure
plot(FourRef, Fref)




sRate=256;    
[FourRef,Fref]=fftBasic(ns_avg,round(sRate));
figure
plot(FourRef, Fref)

pxx2 = pwelch(ns_avg') %richtig
figure
plot(pxx2)


dat=squeeze(mean(avg_data_appended.trial));


% 
% dat2=reshape(dat,248*385,1)
% dat=mean(dat');
% 
pxx = pwelch(dat') %richtig
figure
plot(pxx)

sRate=256;    
[FourRef,Fref]=fftBasic(ns_avg',round(sRate));
figure
plot(FourRef, Fref)









avg_data_appended

sRate=1017;    
[FourRef,Fref]=fftBasic(data_bpfreq.trial{1,1}(1:248, :),round(sRate));
figure
plot(FourRef, Fref)

data_dislike

sRate=256.0001;    
[FourRef,Fref]=fftBasic(vs_allRuns_appended.trial{1,1},round(sRate));
figure
plot(FourRef, Fref)


sRate=256;    
[FourRef,Fref]=fftBasic(data_bpfreq_sel_res.trial{1,1}(1:248, :),round(sRate));
figure
plot(FourRef, Fref)