function [] = adi_prepare_hdm(mniPath, outPath, subj)
if exist ([outPath 'vol.mat'], 'file');
    return
end
if ~exist([mniPath 'mri_realigned.mat'], 'file')
    cfg = [];
    cfg.method = 'interactive';
    cfg.coordsys = 'ctf';
    mri = ft_read_mri([mniPath 'T1warped.mat']);
    mri_realigned = ft_volumerealign(cfg, mri);
    save ([mniPath 'mri_realigned'], 'mri_realigned');
else 
    load ([mniPath 'mri_realigned.mat'])
end

cfg = [];
mri_reslice = ft_volumereslice(cfg, mri_realigned);

%% segment
cfg           = [];
cfg.output    = 'brain';
% cfg.output = {'gray','white','csf','skull','scalp'}; % vermutlich nur für
% EEG relevant
cfg.spmversion = 'spm12';
segmentedmri  = ft_volumesegment(cfg, mri_reslice);
segmentedmri.anatomy = mri_reslice.anatomy;

cfg              = [];
cfg.funparameter = 'brain';
ft_sourceplot(cfg,segmentedmri);


save ([mniPath 'segmentedmri'], 'segmentedmri');


%% head model

cfg = [];
cfg.method = 'singleshell';
vol = ft_prepare_headmodel(cfg, segmentedmri);
save ([outPath 'vol'], 'vol');


end