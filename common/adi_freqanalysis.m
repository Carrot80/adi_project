
%% Hanning taper
cfg              = [];
cfg.output       = 'pow';
cfg.channel      = 'MEG';
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.foi          = 1:1:40;                         % analysis 2 to 30 Hz in steps of 2 Hz 
cfg.t_ftimwin    = ones(length(cfg.foi),1).*1;   % length of time window = 0.5 sec
cfg.toi          = -0.5:0.05:1.5;                  % time window "slides" from -0.5 to 1.5 sec in steps of 0.05 sec (50 ms)
cleanMEG_interp_freq = ft_freqanalysis(cfg, cleanMEG_interp);

cfg = [];
% cfg.baseline = [-0.5 -0.1]; 
% cfg.baselinetype = 'absolute'; 	
% cfg.zlim = [-1.5e-27 1.5e-27];	
cfg.channelname   = 'A177'; % top figure
figure;ft_singleplotTFR(cfg, cleanMEG_interp_freq);

cfg = [];
% cfg.baseline = [-0.5 -0.1]; 
% cfg.baselinetype = 'absolute'; 	
% cfg.zlim = [-1.5e-27 1.5e-27];	
cfg.channelname   = 'A177'; % top figure
figure;ft_singleplotTFR(cfg, cleanMEG_interp_freq);


cfg = [];
cfg.baseline     = [-0.5 -0.1]; 
cfg.baselinetype = 'absolute'; 
cfg.zlim         = [-3e-27 3e-27];	        
cfg.showlabels   = 'yes';	
% cfg.layout       = 'CTF151_helmet.mat';
figure 
ft_multiplotTFR(cfg, cleanMEG_interp_freq);


%% compute sensor level Fourier spectra
cfg            = [];
cfg.method     = 'mtmfft';
cfg.output     = 'powandcsd';
cfg.keeptrials = 'yes';
cfg.baseline = [-0.5 -0.1]; 
cfg.keeptrials = 'yes';
cfg.tapsmofrq  = 2;
% cfg.foi        = foi;
cfg.foilim = [1 40]
cfg.channel     = {'MEG'};
% cfg.trials = cleanMEG_interp.trial;
cfg.taper =  'dpss';
freq           = ft_freqanalysis(cfg, cleanMEG_interp);





cfg = [];
% cfg.baseline = [-0.5 -0.1]; 
% cfg.baselinetype = 'absolute'; 	
% cfg.zlim = [-1.5e-27 1.5e-27];	
% cfg.channelname   = 'MLC'; % top figure
figure;ft_singleplotTFR(cfg, freq);



% cfg = [];
% cfg.method    = 'mtmfft';
% cfg.output    = 'powandcsd'; 
% cfg.tapsmofrq = 2;
% cfg.foilim    = [2 2];
% freqAll = ft_freqanalysis(cfg, dataAll);
% Then we compute the inverse filter based on both conditions. Note the use of cfg.keepfilter so that the output saves this computed filter.
