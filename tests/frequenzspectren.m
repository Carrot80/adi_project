ns_avg=vs_allRuns_appended.trial{1,1};
ns_avg=ns_avg';

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