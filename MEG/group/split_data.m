function [varargout] = adi_split_data(varargin)

if size(varargin,2) > 2
    
    [like_training, dislike_training] = adi_adddata(varargin);
     out.like_training = like_training;
     out.dislike_training = dislike_training;
     
elseif 1 == strcmp(varargin{2}, 'like_dislike')
    
    session = varargin{1};
    clear varargin
    
    ind_like = find(session.labels == 1);
    like = [];
    like.trial = session.trial(ind_like);
    like.time = session.time(ind_like);
    like.response_label = session.response_label(ind_like);
    like.balldesign_short = session.balldesign_short(ind_like);
    like.subject = session.subject(ind_like);
    like.trial_info = session.trial_info(ind_like);
    like.label = session.label;
    like.fsample = session.fsample;
    like.grad = session.grad;
    for k = 1:length(like.subject)
       like.subject_num(k) =  str2double(like.subject{k}(end-1:end));
    end
    
     %% condition dislike:

    ind_dislike = find(session.labels == 2);
    dislike = [];
    dislike.trial = session.trial(ind_dislike);
    dislike.time = session.time(ind_dislike);
    dislike.response_label = session.response_label(ind_dislike);
    dislike.balldesign_short = session.balldesign_short(ind_dislike);
    dislike.trial_info = session.trial_info(ind_dislike);
    dislike.subject = session.subject(ind_dislike);
    dislike.label = session.label;
    dislike.fsample = session.fsample;
    dislike.grad = session.grad;
    for k = 1:length(dislike.subject)
       dislike.subject_num(k) =  str2double(dislike.subject{k}(end-1:end));
    end
    
    varargout{1} = like;
    varargout{2} = dislike;
    
else
    
    session = varargin{1};
    clear varargin
    
    %% condition like:
    % sortiere Probanden nach Trialzahl, absteigend:
    
    ind_like = find(session.labels == 1);
    like = [];
    like.trial = session.trial(ind_like);
    like.time = session.time(ind_like);
    like.response_label = session.response_label(ind_like);
    like.balldesign_short = session.balldesign_short(ind_like);
    like.subject = session.subject(ind_like);
    for k = 1:length(like.subject)
       like.subject_num(k) =  str2double(like.subject{k}(end-1:end));
    end

    [count_trl_like_subject, subj_like] = histcounts(categorical(like.subject_num), categorical(unique(like.subject_num)));
    num_subj = unique(like.subject);
    subj_num_max = find(count_trl_like_subject==max(count_trl_like_subject));

    [num_trl_sorted, ind_sorted]  = sort(count_trl_like_subject,'descend');
    subj_sorted = subj_like(ind_sorted);

    subj_ind_sorted = find(like.subject_num==str2double(subj_sorted(1)));
    for k=2:length(subj_sorted)
        ind = find(like.subject_num==str2double(subj_sorted(k)));
        subj_ind_sorted = cat(2, subj_ind_sorted, ind);
        ind = [];
    end

    like.trial = like.trial(subj_ind_sorted);
    like.time = like.time(subj_ind_sorted);
    like.response_label = like.response_label(subj_ind_sorted);
    like.balldesign_short = like.balldesign_short(subj_ind_sorted);
    like.subject = like.subject(subj_ind_sorted);
    like.subject_num = like.subject_num(subj_ind_sorted);
    order_subj_like_final = unique(like.subject, 'stable');

    %% condition dislike:
    % sortiere Probanden nach Trialzahl, absteigend:
    ind_dislike = find(session.labels == 2);
    dislike = [];
    dislike.trial = session.trial(ind_dislike);
    dislike.time = session.time(ind_dislike);
    dislike.response_label = session.response_label(ind_dislike);
    dislike.balldesign_short = session.balldesign_short(ind_dislike);
    dislike.subject = session.subject(ind_dislike);
    for k = 1:length(dislike.subject)
       dislike.subject_num(k) =  str2double(dislike.subject{k}(end-1:end));
    end

    [count_trl_dislike_subject, subj_dislike] = histcounts(categorical(dislike.subject_num), categorical(unique(dislike.subject_num)));
    num_subj = unique(dislike.subject);
    subj_num_max = find(count_trl_dislike_subject==max(count_trl_dislike_subject));

    [num_trl_sorted, ind_sorted]  = sort(count_trl_dislike_subject,'descend');
    subj_sorted = subj_dislike(ind_sorted);

    subj_ind_sorted = find(dislike.subject_num==str2double(subj_sorted(1)));
    for k=2:length(subj_sorted)
        ind = find(dislike.subject_num==str2double(subj_sorted(k)));
        subj_ind_sorted = cat(2, subj_ind_sorted, ind);
        ind = [];
    end

    dislike.trial = dislike.trial(subj_ind_sorted);
    dislike.time = dislike.time(subj_ind_sorted);
    dislike.response_label = dislike.response_label(subj_ind_sorted);
    dislike.balldesign_short = dislike.balldesign_short(subj_ind_sorted);
    dislike.subject = dislike.subject(subj_ind_sorted);
    dislike.subject_num = dislike.subject_num(subj_ind_sorted);

    order_subj_dislike_final = unique(dislike.subject, 'stable');


    %% cross validation like:
    like_cv = [];
    like_cv.trial = like.trial([find(strcmp(like.subject, char(order_subj_like_final(1)))) find(strcmp(like.subject, char(order_subj_like_final(2))))]);
    like_cv.time = like.time([find(strcmp(like.subject, char(order_subj_like_final(1)))) find(strcmp(like.subject, char(order_subj_like_final(2))))]);
    like_cv.subject = like.subject([find(strcmp(like.subject, char(order_subj_like_final(1)))) find(strcmp(like.subject, char(order_subj_like_final(2))))]);
    like_cv.response_label = like.response_label([find(strcmp(like.subject, char(order_subj_like_final(1)))) find(strcmp(like.subject, char(order_subj_like_final(2))))]);
    like_cv.grad = session.grad;
    like_cv.label = session.label;
    like_cv.labels_cond = ones(1, size(like_cv.trial,2));

    % choose 2 subjects for training
    like_training = [];
    like_training.trial = like.trial([find(strcmp(like.subject, char(order_subj_like_final(3)))) find(strcmp(like.subject, char(order_subj_like_final(4))))]);
    like_training.time = like.time([find(strcmp(like.subject, char(order_subj_like_final(3)))) find(strcmp(like.subject, char(order_subj_like_final(4))))]);
    like_training.subject = like.subject([find(strcmp(like.subject, char(order_subj_like_final(3)))) find(strcmp(like.subject, char(order_subj_like_final(4))))]);
    like_training.response_label = like.response_label([find(strcmp(like.subject, char(order_subj_like_final(3)))) find(strcmp(like.subject, char(order_subj_like_final(4))))]);
    like_training.balldesign_short = like.balldesign_short([find(strcmp(like.subject, char(order_subj_like_final(3)))) find(strcmp(like.subject, char(order_subj_like_final(4))))]);
    like_training.label = session.label;
    like_training.grad = session.grad;

    %% cross validation dislike:

    dislike_cv = [];
    dislike_cv.trial = dislike.trial([find(strcmp(dislike.subject, char(order_subj_dislike_final(1)))) find(strcmp(dislike.subject, char(order_subj_dislike_final(2))))]);
    dislike_cv.time = dislike.time([find(strcmp(dislike.subject, char(order_subj_dislike_final(1)))) find(strcmp(dislike.subject, char(order_subj_dislike_final(2))))]);
    dislike_cv.subject = dislike.subject([find(strcmp(dislike.subject, char(order_subj_dislike_final(1)))) find(strcmp(dislike.subject, char(order_subj_dislike_final(2))))]);
    dislike_cv.response_label = dislike.response_label([find(strcmp(dislike.subject, char(order_subj_dislike_final(1)))) find(strcmp(dislike.subject, char(order_subj_dislike_final(2))))]);
    dislike_cv.grad = session.grad;
    dislike_cv.label = session.label;
    dislike_cv.labels_cond = 2*ones(1, size(like_cv.trial,2));

    % choose 2 subjects for training

    dislike_training = [];
    dislike_training.trial = dislike.trial([find(strcmp(dislike.subject, char(order_subj_dislike_final(3)))) find(strcmp(dislike.subject, char(order_subj_dislike_final(4))))]);
    dislike_training.time = dislike.time([find(strcmp(dislike.subject, char(order_subj_dislike_final(3)))) find(strcmp(dislike.subject, char(order_subj_dislike_final(4))))]);
    dislike_training.subject = dislike.subject([find(strcmp(dislike.subject, char(order_subj_dislike_final(3)))) find(strcmp(dislike.subject, char(order_subj_dislike_final(4))))]);
    dislike_training.response_label = dislike.response_label([find(strcmp(dislike.subject, char(order_subj_dislike_final(3)))) find(strcmp(dislike.subject, char(order_subj_dislike_final(4))))]);
    dislike_training.balldesign_short = dislike.balldesign_short([find(strcmp(dislike.subject, char(order_subj_dislike_final(3)))) find(strcmp(dislike.subject, char(order_subj_dislike_final(4))))]);
    dislike_training.label = session.label;
    dislike_training.grad = session.grad;

%     if length(like_cv.trial) ~= length(dislike_cv.trial)
%         error('unequal number of trials in cross validation set')
%     elseif length(like_training.trial) ~= length(dislike_training.trial)
%          error('unequal number of trials in training set')
%     end

    order_subj_like_final(1:4)=[];
    order_subj_dislike_final(1:4)=[];
    
    out.like_training = like_training;
    out.dislike_training = dislike_training;
    out.like_cv = like_cv;
    out.dislike_cv = dislike_cv;
    out.order_subj_final.like = order_subj_like_final;
    out.order_subj_final.dislike = order_subj_dislike_final; 
    
end
end

function   [like_training, dislike_training] = adi_adddata(vargin)

session = vargin{1};
like_training = vargin{2};
dislike_training = vargin{3};
subj_like = vargin{4};
subj_dislike = vargin{5};
clear vargin
 
 %% condition like:
if ~isempty(subj_like)   
 ind_subj_like = find(strcmp(session.subject, subj_like));
 
 like_training.trial = cat(2,like_training.trial, session.trial(ind_subj_like));
 like_training.time = cat(2,like_training.time, session.time(ind_subj_like));
 like_training.response_label = cat(2,like_training.response_label, session.response_label(ind_subj_like));
 like_training.balldesign_short = cat(2,like_training.balldesign_short, session.balldesign_short(ind_subj_like));
 like_training.subject = cat(2,like_training.subject, session.subject(ind_subj_like));
end
    %% condition dislike:
if ~isempty(subj_dislike)    
 ind_subj_dislike = find(strcmp(session.subject, subj_dislike));
 dislike_training.trial = cat(2,dislike_training.trial, session.trial(ind_subj_dislike));
 dislike_training.time = cat(2,dislike_training.time, session.time(ind_subj_dislike));
 dislike_training.response_label = cat(2,dislike_training.response_label, session.response_label(ind_subj_dislike));
 dislike_training.balldesign_short = cat(2,dislike_training.balldesign_short, session.balldesign_short(ind_subj_dislike));
 dislike_training.subject = cat(2,dislike_training.subject, session.subject(ind_subj_dislike));
end
end