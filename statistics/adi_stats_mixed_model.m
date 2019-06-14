function stats(session, path, balldesign)
% Erstellung der Datenstruktur und ...
% virtuelle Sensoren raussuchen (ohne Kleinhirn, Thalamus, Basalganglien)
load ('U:\My Documents\MATLAB\atlas_clusterstatistic_ohne_cerebellum_thalamus_basalganglien.mat')
% atlas_cluster = sortrows(atlas_cluster, 3);
% for k=1:length(atlas_cluster)
%     atlas_cluster(k,4)={k}; 
% end

ind_cluster = cell2mat(atlas_cluster(:,3));
session.data=session.data(:,ind_cluster, :);

% for k=1:length(like.subject_num)
%     like.subject_num_mat(k)=str2num(like.subject_num{k});
% end
% 
% for k=1:length(dislike.subject_num)
%     dislike.subject_num_mat(k)=str2num(dislike.subject_num{k});
% end

[bincounts,binranges] = histc(str2double(char(session.subject_num)), unique(str2double(char(session.subject_num))));

figure
bar(unique(str2double(char(session.subject_num))),bincounts,'histc')

roi_num = unique(cell2mat(atlas_cluster(:,1)));

for k=1:length(roi_num)
    ind_rois = find(cell2mat(atlas_cluster(:,1))==roi_num(k));
end

for k=1:length(roi_num)
    roi(k).ind = find(cell2mat(atlas_cluster(:,1))==roi_num(k));
    label(k,1) = atlas_cluster(roi(k).ind(1),2);
end

for i = 1:size(session.data,1)
    for k=1:length(roi_num)
        data_roi(i,k,:) = rms(session.data(i, roi(k).ind, :)); 
    end
end

% comp 1: 50-160ms; 
% comp 2: 160-500ms
% comp 3: 500-800ms

time_comp1 = [0.05 0.16];
time_comp2 = [0.16 0.5];
time_comp3 = [0.5 0.8];

sample_comp1 = [nearest(session.time{1,1},time_comp1(1)) nearest(session.time{1,1},time_comp1(2))];
sample_comp2 = [nearest(session.time{1,1},time_comp2(1)) nearest(session.time{1,1},time_comp2(2))];
sample_comp3 = [nearest(session.time{1,1},time_comp3(1)) nearest(session.time{1,1},time_comp3(2))];

for i = 1:length(data_roi)
    data_rois_comp1(i,:) = squeeze(mean(data_roi(i,:, sample_comp1(1):sample_comp1(2)),3));
    data_rois_comp2(i,:) = squeeze(mean(data_roi(i,:, sample_comp2(1):sample_comp2(2)),3));
    data_rois_comp3(i,:) = squeeze(mean(data_roi(i,:, sample_comp3(1):sample_comp3(2)),3));
end

data.subject = session.subject_name';
data.balldesign = session.balldesign';

for i = 1:size(label,1)
    data.([label{i} '_comp1']) = data_rois_comp1(:,i);
    data.([label{i} '_comp2']) = data_rois_comp2(:,i);
    data.([label{i} '_comp3']) = data_rois_comp3(:,i);
end

data.response = session.labels';

data_table = struct2table(data);

adi_glm (data_table);


end

function [] = plot_conditions (session)
ind_like=find(session.labels==1);
ind_dislike=find(session.labels==2);

% like.data = session.data(ind_like,ind_cluster,:);
like.data = session.data(ind_like,:,:);
like.subject_num = session.subject_num(ind_like);
% dislike.data = session.data(ind_dislike,ind_cluster,:);
dislike.subject_num =  session.subject_num(ind_dislike);
dislike.data = session.data(ind_dislike,:,:);

like_rms = mean(squeeze(rms(like.data)));
dislike_rms = mean(squeeze(rms(dislike.data)));
figure
plot(session.time{1,1}, like_rms)
hold on
plot(session.time{1,1},dislike_rms)
legend({'like'; 'dislike'})


end


function [] = alt ()



% L=length(roi_num).*size(session.data,1);
% data_rois_comp1_reshape = reshape(data_rois_comp1, L,1);
% data_rois_comp2_reshape = reshape(data_rois_comp2, L,1);
% data_rois_comp3_reshape = reshape(data_rois_comp3, L,1);

% roi_labels = repmat(label, size(session.data,1),1);
% roi_num_rep = repmat(roi_num, size(session.data,1),1);

% rep_subj=repmat(session.subject_name(1), length(roi_num),1);
% for p=2:length(session.subject_name)
%     rep_subj_s=repmat(session.subject_name(p), length(roi_num),1);
%     rep_subj = cat(1, rep_subj, rep_subj_s);
% end


rep_subj_ID=repmat(str2double(session.subject_num{1}), length(roi_num),1);
for p=2:length(session.subject_num)
    rep_subj_ID_s=repmat(str2double(session.subject_num{p}), length(roi_num),1);
    rep_subj_ID = cat(1, rep_subj_ID, rep_subj_ID_s);
end

rep_balldesign=repmat(session.balldesign(1), length(roi_num),1);
for p=2:length(session.balldesign)
    rep_balldesign_s=repmat(session.balldesign{p}, length(roi_num),1);
    rep_balldesign = cat(1, rep_balldesign, rep_balldesign_s);
end

rep_response=repmat(session.labels(1), length(roi_num),1);
for p=2:length(session.labels)
    rep_response_s=repmat(session.labels(p), length(roi_num),1);
    rep_response = cat(1, rep_response, rep_response_s);
end

for p = 1:L
    data(p).subject = rep_subj(p);
    data(p).ID = rep_subj_ID(p);
    data(p).balldesign = rep_balldesign{p};
    data(p).roi= roi_labels{p};
    data(p).comp1 = data_rois_comp1_reshape(p);
    data(p).comp2 = data_rois_comp2_reshape(p);
    data(p).comp3 = data_rois_comp3_reshape(p);
    data(p).response = rep_response(p);
end







end