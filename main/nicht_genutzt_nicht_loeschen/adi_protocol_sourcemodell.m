

%% Prepare mri and headmodel
if ~exist([parentPathOut 'vol.mat'], 'file')
    [vol, mri, segmentedmri] = conn_prepare_hdm(mripath, headshapepath);
    save([parentPathOut 'vol.mat'], 'vol');
    save([parentPathOut 'mri.mat'], 'mri');
    save([parentPathOut 'segmentedmri.mat'], 'segmentedmri');
    disp('Completed MRI preparation.');
else
    disp('Files exist, skipping processing step.');
    load([parentPathOut 'vol.mat']);
    load([parentPathOut 'mri.mat']);
    load([parentPathOut 'segmentedmri.mat']);
end

%% Prepare sourcemodel
savepath = fullfile(parentPathOut, 'sourcemodel_template.mat');
if ~exist(savepath, 'file')
    sourcemodel = conn_prepare_sourcemodel(mri, vol, 1, 0);
    save(savepath, 'sourcemodel');
    disp('Completed sourcemodel.');
else
    disp('Files exist, skipping processing step.');
    load(savepath);
end

%% Prepare data
savepath = fullfile(outpath, ['dataclean' num2str(triggercode) '.mat']);
if ~exist(savepath, 'file')
    dataclean = conn_prepare_data_trials(datapath, triggercode, pretrigger, posttrigger, megSystem);
    save(savepath, 'dataclean', '-v7.3');
    disp('Completed data preparation.');
else
    disp('Files exist, skipping processing step.');
    load([outpath 'dataclean.mat']);
end

%% Denoise
savepath = fullfile(outpath, ['datadenoised' num2str(triggercode) '.mat']);

if ~exist(savepath, 'file')
%     cfg = [];
%     cfg.bpfilter = 'yes';
%     cfg.bpfreq = [4 7];
%     filtered = ft_preprocessing(cfg, dataclean);
    
    values = nt_trial2mat(dataclean.trial);
    c0 = nt_cov(values);
    c1 = nt_cov(mean(values, 3));
    todss = nt_dss0(c0, c1);
    z = nt_mmat(values, todss);
    clean = nt_tsr(values, z(:,8:end, :));
    denoised = dataclean;
    denoised.trial = nt_mat2trial(clean);
    
    save(savepath, 'denoised', '-v7.3');
    disp('Completed data denoising.');
else
    disp('Files exist, skipping processing step.');
    load(savepath);
end

%% Select data
cfg = [];
cfg.toilim = [-0.8 -0.2];
data_pre = ft_redefinetrial(cfg, denoised);
cfg.toilim = [0.2 0.8];
data_post = ft_redefinetrial(cfg, denoised);

%% Analyze
filepath = fullfile(outpath, 'results_delta_pre.mat');
if ~exist(filepath, 'file')
    [results_delta_pre, source_conn_delta_pre] = conn_analyze_freqs(data_pre, vol, sourcemodel, 1:3);
    save(filepath, 'results_delta_pre', 'source_conn_delta_pre');
    disp('Completed delta pre analysis.');
else
    disp('Files exist, skipping processing step.');
    load(filepath);
end

filepath = fullfile(outpath, 'results_delta_post.mat');
if ~exist(filepath, 'file')
    [results_delta_post, source_conn_delta_post] = conn_analyze_freqs(data_post, vol, sourcemodel, 1:3);
    save(filepath, 'results_delta_post', 'source_conn_delta_post');
    disp('Completed delta post analysis.');
else
    disp('Files exist, skipping processing step.');
    load(filepath);
end

filepath = fullfile(outpath, 'results_theta_pre.mat');
if ~exist(filepath, 'file')
    [results_theta_pre, source_conn_theta_pre] = conn_analyze_freqs(data_pre, vol, sourcemodel, 4:7);
    save(filepath, 'results_theta_pre', 'source_conn_theta_pre');
    disp('Completed theta pre analysis.');
else
    disp('Files exist, skipping processing step.');
    load(filepath);
end

filepath = fullfile(outpath, 'results_theta_post.mat');
if ~exist(filepath, 'file')
    [results_theta_post, source_conn_theta_post] = conn_analyze_freqs(data_post, vol, sourcemodel, 4:7);
    save(filepath, 'results_theta_post', 'source_conn_theta_post');
    disp('Completed theta post analysis.');
else
    disp('Files exist, skipping processing step.');
    load(filepath);
end

filepath = fullfile(outpath, 'results_alpha_pre.mat');
if ~exist(filepath, 'file')
    [results_alpha_pre, source_conn_alpha_pre] = conn_analyze_freqs(data_pre, vol, sourcemodel, 8:15);
    save(filepath, 'results_alpha_pre', 'source_conn_alpha_pre');
    disp('Completed alpha pre analysis.');
else
    disp('Files exist, skipping processing step.');
    load(filepath);
end

filepath = fullfile(outpath, 'results_alpha_post.mat');
if ~exist(filepath, 'file')
    [results_alpha_post, source_conn_alpha_post] = conn_analyze_freqs(data_post, vol, sourcemodel, 8:15);
    save(filepath, 'results_alpha_post', 'source_conn_alpha_post');
    disp('Completed alpha post analysis.');
else
    disp('Files exist, skipping processing step.');
    load(filepath);
end

filepath = fullfile(outpath, 'results_gamma_pre.mat');
if ~exist(filepath, 'file')
    [results_gamma_pre, source_conn_gamma_pre] = conn_analyze_freqs(data_pre, vol, sourcemodel, 25:50);
    save(filepath, 'results_gamma_pre', 'source_conn_gamma_pre');
    disp('Completed gamma pre analysis.');
else
    disp('Files exist, skipping processing step.');
    load(filepath);
end

filepath = fullfile(outpath, 'results_gamma_post.mat');
if ~exist(filepath, 'file')
    [results_gamma_post, source_conn_gamma_post] = conn_analyze_freqs(data_post, vol, sourcemodel, 25:50);
    save(filepath, 'results_gamma_post', 'source_conn_gamma_post');
    disp('Completed gamma post analysis.');
else visual_source_analysis
    disp('Files exist, skipping processing step.');
    load(filepath);
end

filepath = fullfile(outpath, 'results_highgamma_pre.mat');
if ~exist(filepath, 'file')
    [results_highgamma_pre, source_conn_highgamma_pre] = conn_analyze_freqs(data_pre, vol, sourcemodel, 70:100);
    save(filepath, 'results_highgamma_pre', 'source_conn_highgamma_pre');
    disp('Completed highgamma pre analysis.');
else
    disp('Files exist, skipping processing step.');
    load(filepath);
end

filepath = fullfile(outpath, 'results_highgamma_post.mat');
if ~exist(filepath, 'file')
    [results_highgamma_post, source_conn_highgamma_post] = conn_analyze_freqs(data_post, vol, sourcemodel, 70:100);
    save(filepath, 'results_highgamma_post', 'source_conn_highgamma_post');
    disp('Completed highgamma post analysis.');
else
    disp('Files exist, skipping processing step.');
    load(filepath);
end

%%
cfg = [];
cfg.parameter = 'degrees';
cfg.operation = 'x1-x2';
ratio_delta = ft_math(cfg, results_delta_post, results_delta_pre);
ratio_theta = ft_math(cfg, results_theta_post, results_theta_pre);
ratio_alpha = ft_math(cfg, results_alpha_post, results_alpha_pre);
ratio_gamma = ft_math(cfg, results_gamma_post, results_gamma_pre);
ratio_highgamma = ft_math(cfg, results_highgamma_post, results_highgamma_pre);

ratio_theta = results_theta_post;
ratio_theta.degrees = results_theta_post.degrees - results_theta_pre.degrees;

%%
conn_visualize_network(ratio_delta, sourcemodel, mri, 'Delta ratio');
conn_visualize_network(ratio_theta, sourcemodel, mri, 'Theta ratio');
conn_visualize_network(ratio_alpha, sourcemodel, mri, 'Alpha ratio');
conn_visualize_network(ratio_gamma, sourcemodel, mri, 'Gamma ratio');
conn_visualize_network(ratio_highgamma, sourcemodel, mri, 'High gamma ratio');

%% Visualize results
conn_visualize_network(results_delta_pre, sourcemodel, mri, 'Delta pre');
conn_visualize_network(results_delta_post, sourcemodel, mri, 'Delta post');
results_delta_diff = results_delta_post;
results_delta_diff.degrees = results_delta_diff.degrees - results_delta_pre.degrees;
conn_visualize_network(results_delta_diff, sourcemodel, mri, 'Delta diff');

conn_visualize_network(results_theta_pre, sourcemodel, mri, 'Theta pre');
conn_visualize_network(results_theta_post, sourcemodel, mri, 'Theta post');
results_theta_diff = results_theta_post;
results_theta_diff.degrees = results_theta_diff.degrees - results_theta_pre.degrees;
conn_visualize_network(results_theta_diff, sourcemodel, mri, 'Theta diff');

conn_visualize_network(results_alpha_pre, sourcemodel, mri, 'Alpha pre');
conn_visualize_network(results_alpha_post, sourcemodel, mri, 'Alpha post');
results_alpha_diff = results_alpha_post;
results_alpha_diff.degrees = results_alpha_diff.degrees - results_alpha_pre.degrees;
conn_visualize_network(results_alpha_diff, sourcemodel, mri, 'Alpha diff');

conn_visualize_network(results_gamma_pre, sourcemodel, mri, 'Gamma pre');
conn_visualize_network(results_gamma_post, sourcemodel, mri, 'Gamma post');
results_gamma_diff = results_gamma_post;
results_gamma_diff.degrees = results_gamma_diff.degrees - results_gamma_pre.degrees;
conn_visualize_network(results_gamma_diff, sourcemodel, mri, 'Gamma diff');

results_highgamma_diff = results_highgamma_post;
results_highgamma_diff.degrees = results_highgamma_diff.degrees - results_highgamma_pre.degrees;
conn_visualize_network(results_highgamma_diff, sourcemodel, mri, 'High Gamma diff');

%% Parcelate
parcel_theta = conn_parcellate(ratio_theta, 'degrees');
parcel_alpha = conn_parcellate(ratio_alpha, 'degrees');


%% Laterality
lat_theta = conn_roiact(parcel_theta, 'degrees');
lat_alpha = conn_roiact(parcel_alpha, 'degrees');
