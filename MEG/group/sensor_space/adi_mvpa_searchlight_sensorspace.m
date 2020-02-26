function stats(like, dislike, time, path2save, ball)
%% searchlight analyse: where

% comp: Zeitintervall extrahiert aus rms

cfg = [];
cfg.method      = 'triangulation';
cfg.layout      = '4D248.lay';    
cfg.channel     = dislike.label;
neighbours = ft_prepare_neighbours(cfg);



%%
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
    title(['comp ' num2str(time(1)) ' ' num2str(time(2)) ' ' ball]);
    
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
