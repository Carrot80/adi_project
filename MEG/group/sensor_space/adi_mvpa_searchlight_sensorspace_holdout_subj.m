function stats(like, dislike, time, path2save, ball)
%% searchlight analyse: where
rng_state = rng;

cfg = [];
cfg.method      = 'triangulation';
cfg.layout      = '4D248.lay';    
cfg.channel     = dislike.label;
neighbours = ft_prepare_neighbours(cfg);
neighbours = adi_channelconnectivity(struct('neighbours',neighbours, 'channel', {cfg.channel}));

% figure,
% imagesc(neighbours)
% set(gca,'XTickLabel',like.label(get(gca,'XTick')))
% set(gca,'YTickLabel',like.label(get(gca,'YTick')))
% title('Neighbourhood matrix')
% grid on

% comp: Zeitintervall extrahiert aus rms
t1_sample = nearest(like.time{1}, time(1));
t2_sample = nearest(like.time{1}, time(2));
like.trial = like.trial(:,:, t1_sample : t2_sample);
dislike.trial = dislike.trial(:,:, t1_sample : t2_sample);

like.trial = mean(like.trial(:,:, :),3);
dislike.trial = mean(dislike.trial(:,:, :),3);

CV = adi_crossval(like.subject_num, dislike.subject_num, 'holdout');
mvpa = {'lda'; 'logreg'; 'svm'}; 
   
    %% built train and test data based on CV folds
cf = cell(length(mvpa),length(CV.pairs));
    
for cc = 1:length(CV.pairs)
    
    if 2 == length(CV.like(CV.pairs(cc,1)).test)  && 1 == length(CV.dislike(CV.pairs(cc,2)).test)
        X_test  = cat(1, like.trial(find(ismember(like.subject_num, CV.like(CV.pairs(cc,1)).test(1))),:), ...
                  like.trial(find(ismember(like.subject_num, CV.like(CV.pairs(cc,1)).test(2))),:), ...
                  dislike.trial(find(ismember(dislike.subject_num, CV.dislike(CV.pairs(cc,2)).test(1))),:));    
        
        clabel_test = [ones(1, size(like.trial(find(ismember(like.subject_num, CV.like(CV.pairs(cc,1)).test(1)))),2)) ...
                  ones(1, size(like.trial(find(ismember(like.subject_num, CV.like(CV.pairs(cc,1)).test(2)))),2)) ...
                  2*ones(1,size(dislike.trial(find(ismember(dislike.subject_num, CV.dislike(CV.pairs(cc,2)).test(1)))),2))];
        X_train_like = like.trial;
        X_train_like([find(ismember(like.subject_num, CV.like(CV.pairs(cc,1)).test(1))) find(ismember(like.subject_num, CV.like(CV.pairs(cc,1)).test(2)))],:) = []; 
        X_train_dislike = dislike.trial;
        X_train_dislike(find(ismember(dislike.subject_num, CV.dislike(CV.pairs(cc,2)).test(1))),:) = []; 

    elseif 2 == length(CV.dislike(CV.pairs(cc,2)).test) && 1 == length(CV.like(CV.pairs(cc,1)).test)
        X_test  = cat(1, like(find(ismember(like.subject_num, CV.like(CV.pairs(cc,1)).test(1))),:), ...
                  dislike.trial(find(ismember(dislike.subject_num, CV.dislike(CV.pairs(cc,2)).test(1))),:), ...
                  dislike.trial(find(ismember(dislike.subject_num, CV.dislike(CV.pairs(cc,2)).test(2))),:));
        
        clabel_test = [ones(1, size(like(find(ismember(like.subject_num, CV.like(CV.pairs(cc,1)).test(1)))),2)) ...
                      2*ones(1,size(dislike.trial(find(ismember(dislike.subject_num, CV.dislike(CV.pairs(cc,2)).test(1)))),2))...
                      2*ones(1,size(dislike.trial(find(ismember(dislike.subject_num, CV.dislike(CV.pairs(cc,2)).test(2)))),2))];
                  
        X_train_like = like.trial;
        X_train_like(find(ismember(like.subject_num, CV.like(CV.pairs(cc,1)).test(1))),:) = []; 
        X_train_dislike = dislike.trial;
        X_train_dislike([find(ismember(dislike.subject_num, CV.dislike(CV.pairs(cc,2)).test(1))) find(ismember(dislike.subject_num, CV.dislike(CV.pairs(cc,2)).test(2)))],:) = []; 

    elseif 1 == length(CV.dislike(CV.pairs(cc,2)).test) && 1 == length(CV.like(CV.pairs(cc,1)).test)
        
        X_test  = cat(1, like(find(ismember(like.subject_num, CV.like(CV.pairs(cc,1)).test(1))),:), ...
                  dislike.trial(find(ismember(dislike.subject_num, CV.dislike(CV.pairs(cc,2)).test(1))),:));
         
        clabel_test = [ones(1, size(like.trial(find(ismember(like.subject_num, CV.like(CV.pairs(cc,1)).test(1)))),2)) ... ...
                  2*ones(1,size(dislike.trial(find(ismember(dislike.subject_num, CV.dislike(CV.pairs(cc,2)).test(1)))),2))];
 
        X_train_like = like.trial;
        X_train_like(find(ismember(like.subject_num, CV.like(CV.pairs(cc,1)).test(1))),:) = []; 
        X_train_dislike = dislike.trial;
        X_train_dislike(find(ismember(dislike.subject_num, CV.dislike(CV.pairs(cc,2)).test(1))),:) = []; 
   
    else
        X_test  = cat(1, like.trial(find(ismember(like.subject_num, CV.like(CV.pairs(cc,1)).test(1))),:), ...
                  like.trial(find(ismember(like.subject_num, CV.like(CV.pairs(cc,1)).test(2))),:), ...
                  dislike.trial(find(ismember(dislike.subject_num, CV.dislike(CV.pairs(cc,2)).test(1))),:), ...
                  dislike.trial(find(ismember(dislike.subject_num, CV.dislike(CV.pairs(cc,2)).test(2))),:));
        
        clabel_test = [ones(1, size(like.trial(find(ismember(like.subject_num, CV.like(CV.pairs(cc,1)).test(1)))),2)) ...
                  ones(1, size(like.trial(find(ismember(like.subject_num, CV.like(CV.pairs(cc,1)).test(2)))),2)) ...
                  2*ones(1,size(dislike.trial(find(ismember(dislike.subject_num, CV.dislike(CV.pairs(cc,2)).test(1)))),2))...
                  2*ones(1,size(dislike.trial(find(ismember(dislike.subject_num, CV.dislike(CV.pairs(cc,2)).test(2)))),2))];
        
        X_train_like = like.trial;
        X_train_like([find(ismember(like.subject_num, CV.like(CV.pairs(cc,1)).test(1))) find(ismember(like.subject_num, CV.like(CV.pairs(cc,1)).test(2)))],:) = []; 
        X_train_dislike = dislike.trial;
        X_train_dislike([find(ismember(dislike.subject_num, CV.dislike(CV.pairs(cc,2)).test(1))) find(ismember(dislike.subject_num, CV.dislike(CV.pairs(cc,2)).test(2)))],:) = []; 
   
    end

    X_train = cat(1, X_train_like, X_train_dislike) ;
    clabel_train = [ones(1, size(X_train_like,1)) 2*ones(1, size(X_train_dislike,1))];
    clear X_train_like X_train_dislike   

    %%  balance data    
    [X_train, clabel_train, labelidx_train] = mv_balance_classes(X_train,clabel_train,'undersample',[]);
    [X_test, clabel_test, labelidx_test] = mv_balance_classes(X_test,clabel_test,'undersample',[]);
    
    for kk = 1 : length(mvpa)
        train_fun = eval(['@train_' mvpa{kk}]);
        test_fun = eval(['@test_' mvpa{kk}]);
        param.(mvpa{kk}) =  mv_get_classifier_param(mvpa{kk});
        
        cf{kk, cc} = train_fun(param.(mvpa{kk}), X_train, clabel_train);
        [predlabel, dval, prob] = train_fun(cf{kk, cc}, X_test);
        % To calculate classification accuracy, compare the predicted labels to
        % the true labels and take the mean
        
        % Calculate AUC and Accuracy
        perf.(mvpa{kk}).auc(cc) = mv_calculate_performance('auc', 'dval', dval, clabel_test);
        perf.(mvpa{kk}).accuracy(cc) = mv_calculate_performance('acc', 'dval', dval, clabel_test);

    end
cfg = [];
cfg.method      = 'triangulation';
cfg.layout      = '4D248.lay';    
cfg.channel     = dislike.label;
neighbours = ft_prepare_neighbours(cfg);
rng default
time_ = [num2str(time(1)) '_' num2str(time(2)) 's'];
mvpa = {'lda'; 'logreg'}; 

for p = 1:length(mvpa)

    cfg = [] ;
    cfg.method          = 'mvpa'; 
    cfg.searchlight = 'yes';
    cfg.mvpa.neighbours  = neighbours;
    cfg.mvpa.classifier = mvpa{p}; %  % logreg, lda, svm, ensemble, kernel_fda (for nonlinear data), libsvm, liblinear  multi-class Linear Discriminant Analysis (LDA)
    cfg.mvpa.metric     = {'accuracy'; 'auc'};
    % cfg.mvpa.param.lambda =  [0.0001 0.0006 0.0036 0.026 0.13 0.1 0.2 0.3 0.4 0.5 0.6 0.77 1 4.6416 27.82559 166.81 1000];
    % cfg.mvpa.param.reg = 'l2';%'l2';
    cfg.mvpa.param.plot = 0;
    cfg.mvpa.param.repeat = 5;
    cfg.mvpa.balance = 'undersample';
    cfg.mvpa.repeat = 5;
    cfg.mvpa.normalise = 'none';
    cfg.latency     = time; 
    cfg.avgovertime = 'yes';
    cfg.design          = [ones(1, size(like.trial,2)) 2*ones(1,size(dislike.trial,2))]';
    stat_comp.(mvpa{p}) = ft_timelockstatistics(cfg, like, dislike);
    save([path2save 'stat_comp_' time_ '.mat'], 'stat_comp'); 
%     fprintf('Classification accuracy: %0.2f\n', stat_comp.accuracy)

    figure
    cfg_fig              = [];
    cfg_fig.parameter    = 'accuracy'; %oder accuracy
    cfg_fig.layout       = '4D248.lay';            
    cfg_fig.xlim         = [0, 0];
    cfg_fig.colorbar     = 'yes';
    ft_topoplotER(cfg_fig, stat_comp.(mvpa{p}));
    title(['comp ' num2str(time(1)) ' ' num2str(time(2))]);
    
    savefig([path2save 'searchlight_accuracy_' mvpa{p} '_' time_ '.fig'])
    close
    
    figure
    cfg_fig              = [];
    cfg_fig.parameter    = 'auc'; 
    cfg_fig.layout       = '4D248.lay';            
    cfg_fig.xlim         = [0, 0];
    cfg_fig.colorbar     = 'yes';
    ft_topoplotER(cfg_fig, stat_comp.(mvpa{p}));
    title(['comp ' num2str(time(1)) ' ' num2str(time(2))]);
    
    savefig([path2save 'searchlight_auc_' mvpa{p} '_' time_ '.fig'])
    close
    
end

end
