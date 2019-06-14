load 'U:\My Documents\adidas_Kirsten\Tutorials\mvpa\dataFIC_LP.mat';
load 'U:\My Documents\adidas_Kirsten\Tutorials\mvpa\dataFC_LP.mat';
load 'U:\My Documents\adidas_Kirsten\Tutorials\mvpa\dataIC_LP.mat';

nFIC = numel(dataFIC_LP.trial);
nFC = numel(dataFC_LP.trial);
nIC = numel(dataIC_LP.trial);

cfg = [] ;
cfg.method          = 'mvpa';
cfg.mvpa.classifier = 'logreg'; %  multi-class Linear Discriminant Analysis (LDA)
cfg.mvpa.metric     = 'auc';
cfg.mvpa.k          = 3;
cfg.mvpa.param.reg = 'l2';
cfg.mvpa.param.lambda = 'auto';
cfg.mvpa.param.plot =1;
cfg.latency         = [0.5, 0.7];
cfg.avgovertime     = 'yes';
cfg.design          = [ones(nFIC,1); 2*ones(nFC,1);];

stat = ft_timelockstatistics(cfg, dataFIC_LP, dataFC_LP)

dat = ft_appenddata([], dataFIC_LP, dataFC_LP, dataIC_LP);
stat = ft_timelockstatistics(cfg, dat);

fprintf('Classification accuracy: %0.2f\n', stat.accuracy)

cfg.mvpa.metric      = 'confusion';
stat = ft_timelockstatistics(cfg, dataFIC_LP, dataFC_LP, dataIC_LP)

stat.confusion

mv_plot_result(stat.mvpa)

cfg = [] ;  
cfg.method           = 'mvpa';
cfg.mvpa.classifier  = 'lda';
cfg.mvpa.metric      = 'auc';
cfg.mvpa.k           = 10;
cfg.mvpa.repeat      = 2;
cfg.avgovertime = 'no';
cfg.design           = [ones(nFIC,1); 2*ones(nFC,1)];

stat = ft_timelockstatistics(cfg, dataFIC_LP, dataFC_LP);

figure
plot(stat.auc)

mv_plot_result(stat.mvpa, dataFC_LP.time{1})
%% searchlight analysis:
cfg = [] ;  
cfg.method      = 'mvpa';
cfg.searchlight = 'yes';
cfg.mvpa.classifier  = 'logreg';
cfg.design      = [ones(nFIC,1); 2*ones(nFC,1)];
cfg.latency     = [0.3, 0.7];
cfg.avgovertime = 'yes'; % no auch möglich
stat = ft_timelockstatistics(cfg, dataFIC_LP, dataFC_LP);

cfg              = [];
cfg.parameter    = 'accuracy';
cfg.layout       = 'CTF151_helmet.mat';            
cfg.xlim         = [0, 0];
cfg.colorbar     = 'yes';
ft_topoplotER(cfg, stat);shg


cfg = [];
cfg.method      = 'triangulation';
cfg.layout      = 'CTF151_helmet.mat';
cfg.channel     = dataFC_LP.label;
neighbours = ft_prepare_neighbours(cfg);

  cfg = [] ;  
  cfg.method      = 'mvpa';
  cfg.searchlight = 'yes';
  cfg.design      = [ones(nFIC,1); 2*ones(nFC,1)];
  cfg.latency     = [0.3, 0.7];
  cfg.avgovertime = 'yes';

  cfg.mvpa.neighbours  = neighbours;

  stat = ft_timelockstatistics(cfg, dataFIC_LP, dataFC_LP)


  cfg              = [];
  cfg.parameter    = 'accuracy';
  cfg.layout       = 'CTF151_helmet.mat';            
  cfg.xlim         = [0, 0];
  cfg.colorbar     = 'yes';
  ft_topoplotER(cfg, stat);
  
  	cfg = [] ;  
    cfg.method      = 'mvpa';
    cfg.timextime   = 'yes';
    cfg.design      = [ones(nFIC,1); 2*ones(nFC,1)];

    stat = ft_timelockstatistics(cfg, dataFIC_LP, dataFC_LP);

figure
mv_plot_result(stat.mvpa, dataFC_LP.time{1})




