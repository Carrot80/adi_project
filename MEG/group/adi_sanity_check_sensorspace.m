function adi_sanitycheck_sensorspace(session)


cfg=[];
avg = ft_timelockanalysis(cfg, session);
figure
plot(avg.time, avg.avg)

for p=1:length(session.trial)
    ind=find(isnan(session.trial{p}));
    if 1==any(ind)
       error(['trial no ' num2str(p)]) 
    end
    
end

avg=squeeze(mean(session.trial{132}));
sRate=256;    
[FourRef,Fref]=fftBasic(session.trial{13},round(sRate));
figure
plot(FourRef, Fref)

ind=find(isnan(session.data));




like.trial=session.trial(find(session.labels==1));
dislike.trial=session.trial(find(session.labels==2));
like.time=session.time(find(session.labels==1));
dislike.time=session.time(find(session.labels==2));

like.label=session.label;
dislike.label=session.label;




like.label = num2cell([1:length(like.trial{1,1})]);
like2 = like.trial(:);

like.avg = mean(like.trial(:))
like.avg = mean(like.trial(:))




cfg=[];
avg_like = ft_timelockanalysis(cfg, like);
avg_dislike = ft_timelockanalysis(cfg, dislike);
figure
plot(avg_like.time, avg_like.avg)
title('avg_like')
figure
plot(avg_dislike.time, avg_dislike.avg)
title('avg_dislike')
%like
for p = 1:length(like.trial)
    temp = like.trial{1,p};
    like.data(p,:,:) = temp;
    clear temp
end

like.rms_mean=squeeze(mean(rms(like.data)));
figure
plot(avg_like.time, like.rms_mean)

%dislike:
for p = 1:length(dislike.trial)
    temp = dislike.trial{1,p};
    dislike.data(p,:,:) = temp;
    clear temp
end

dislike.rms_mean=squeeze(mean(rms(dislike.data)));
hold on
plot(avg_like.time, dislike.rms_mean)
legend('like', 'dislike')
title('rms like vs dislike')

% 
figure
boxplot([like.rms_mean(1:129), dislike.rms_mean(1:129)],'Notch','on','Labels',{'rms like','rms dislike'},'Whisker',1)
[p,h,stats] = ranksum(like.rms_mean(1:129),dislike.rms_mean(1:129));

figure
boxplot([like.rms_mean(130:end), dislike.rms_mean(130:end)],'Notch','on','Labels',{'rms like','rms dislike'},'Whisker',1)
[p,h,stats] = ranksum(like.rms_mean(130:end),dislike.rms_mean(130:end));
title('poststim-intervall')

time1=nearest(avg_like.time, 0.57);
time2=nearest(avg_like.time, 0.7);
figure
boxplot([like.rms_mean(time1:time2), dislike.rms_mean(time1:time2)],'Notch','on','Labels',{'rms like','rms dislike'},'Whisker',1)
[p,h,stats] = ranksum(like.rms_mean(time1:time2),dislike.rms_mean(time1:time2));
title([num2str(0.57) ' bis ' num2str(0.7) ' s'])






end