function [] = adi_grandavg_timeintervall(path2balldesigns, balldesigns, balldesigns_long)

load ('U:\My Documents\MATLAB\atlas_clusterstatistic_ohne_cerebellum_thalamus_basalganglien.mat')
ind_cluster = cell2mat(atlas_cluster(:,3));

for k=1:length(balldesigns)
    if 0==exist('session', 'var') || 1==exist('session', 'var') && 0==strcmp(session.balldesign{1,1}, balldesigns_long{k}) 
        load ([path2balldesigns balldesigns{k} '_ball\session_' balldesigns{k}])
    end
time=session.time{1,1};    
ind_like=find(session.labels==1);
ind_dislike=find(session.labels==2);
if strcmp(balldesigns{k}, 'gbf') 
    like.data = session.data(ind_like,:,:);
    dislike.data = session.data(ind_dislike,:,:);
else
    like.data = session.data(ind_like,ind_cluster,:);
    dislike.data = session.data(ind_dislike,ind_cluster,:); 
end

rms_.(balldesigns{k}).like = mean(squeeze(rms(like.data)));
rms_.(balldesigns{k}).dislike = mean(squeeze(rms(dislike.data)));
clear session
end

for k=1:length(balldesigns)

    h(k)=figure
    plot(time, rms_.(balldesigns{k}).like)
    hold on
    plot(time,rms_.(balldesigns{k}).dislike)
    h_legend = legend({'like', 'dislike'}, 'Location','Northeast')  ;  
    title(balldesigns{k})
    savefig(h(k), [path2balldesigns balldesigns{k} '_ball\rms_' balldesigns{k} '.fig'])

end
close all

% compute grandavg_rms    
grand_avg_rms_like = (rms_.(balldesigns{1}).like+rms_.(balldesigns{2}).like+rms_.(balldesigns{3}).like+rms_.(balldesigns{4}).like+rms_.(balldesigns{5}).like+rms_.(balldesigns{6}).like+rms_.(balldesigns{7}).like+rms_.(balldesigns{8}).like+rms_.(balldesigns{9}).like)/9;
grand_avg_rms_dislike = (rms_.(balldesigns{1}).dislike+rms_.(balldesigns{2}).dislike+rms_.(balldesigns{3}).dislike+rms_.(balldesigns{4}).dislike+rms_.(balldesigns{5}).dislike+rms_.(balldesigns{6}).dislike+rms_.(balldesigns{7}).dislike+rms_.(balldesigns{8}).dislike+rms_.(balldesigns{9}).dislike)/9;
    
figure
plot(time, grand_avg_rms_like)
hold on
plot(time,grand_avg_rms_dislike)
h_legend = legend({'like', 'dislike'}, 'Location','Northeast')  ;  
title('rms_grand_avg')


    
    






end









