
 function [virtsens_all_subjects]= main(virtsens_all_subjects, path2data, freqname, condition, atlas, i)

 [session] = mk_SVM_struct([], path2data, freqname, condition, atlas);
 
 virtsens_all_subjects(i) = session;
 clear session
 
 end
 
 
 function [session] = mk_SVM_struct(session, path2data, freqname, condition, atlas)
 
 fileName = [path2data '\noROIs\runs_appended\virtsens\virtsens_all_conditions_allRuns_', freqname, '.mat'];
 if ~exist(fileName, 'file')
    [vs_allRuns] = adi_appendvirtsensMEG(path2data, freqname, condition);
 else 
    load (fileName)
 end

 
num_conditions = size(fields(vs_allRuns),1);
fieldnames =   fields(vs_allRuns);       

if 3 == num_conditions
    data = cat(2, vs_allRuns.(condition{1}).trial, vs_allRuns.(condition{2}).trial, vs_allRuns.(condition{3}).trial);
    session.labels = cat(2, ones(1,length(vs_allRuns.(condition{1}).trial)), 2*ones(1, length(vs_allRuns.(condition{2}).trial)), 3*ones(1, length(vs_allRuns.(condition{3}).trial)));
elseif 2 == num_conditions
    data = cat(2, vs_allRuns.(fieldnames{1}).trial, vs_allRuns.(fieldnames{2}).trial);
    switch (fieldnames{1})
        case 'like'
            ind_field_1 = 1;
        case 'dislike' 
            ind_field_1 = 2;
        case 'dontcare'
            ind_field_1 = 3;
    end
    switch (fieldnames{2})
        case 'like'
            ind_field_2 = 1;
        case 'dislike' 
            ind_field_2 = 2;
        case 'dontcare'
            ind_field_2 = 3;
    end       
    session.labels = cat(2, ind_field_1*ones(1, length(vs_allRuns.(fieldnames{1}).trial)), ind_field_2*ones(1, length(vs_allRuns.(fieldnames{2}).trial)));
end
for i = 1:length(data)
    temp = data{1,i};
    session.data(i,:,:) = temp;
end
clear data
session.data = session.data(:,cell2mat(atlas(:,3)),:);
session.time = vs_allRuns.(fieldnames{1}).time;

end
  


     
 %%
function [] = SVM_Guggenmos_tutorial(sessions, outPath, freq, condition, atlas)

% delete cerebellum
atlas_indices = cell2mat(atlas(:,3));
sessions.data = sessions.data(:, atlas_indices,:);

% We set a seed, in order to make analyses reproducible:
rng('default');
rng(10);

% Now we set some parameters. Only the number of permutations and 
% the number of pseudo-trials are free parameters. The number of conditions, 
% sensors, time points and sessions are derived from the data (i.e., from the sessions variable above).

% Parameters
n_perm = 20;  % number of permutations
n_pseudo = 5;  % number of pseudo-trials % evtl. erhöhen auf 10?
n_conditions = length(unique(sessions(1).labels));
n_sensors = size(sessions(1).data, 2);
n_time = size(sessions(1).data, 3);
n_sessions = length(sessions);

% The analytic logic is contained in a nested for loop, with loops for the number of sessions, 
% number of permutations, number of timepoints, number of conditions, and number of conditions again. 
% Overall, the logic contains 4 crucial steps:

% 1. Compute pseudo-trials for the training and test data
% 2. Whiten the training data (here using the Epoch method, which is recommended in our manuscript)
% 3. Fit the classifier to the training data
% 4. Compute classification accuracy on test data
% 

% pre-load mechanism, for convenience
% preload_result = true; % for recomputing the decoding analyses, set to false
% if preload_result
%     load(fullfile(pwd, 'result_decoding.mat'))
% else
    % We define three classifiers that will be compared: Support Vector Machine, Weighted Robust Distance,
    % Gaussian Naive Bayes
    clfs = {'svm', 'gnb', 'weird'};
    for c = 1:length(clfs)
        result.(clfs{c}) = nan(n_sessions, n_perm, n_conditions, n_conditions, n_time);
    end
    for s = 1:n_sessions

        fprintf('sessions %g / %g\n', s, n_sessions)

        X = sessions(s).data;
        y = sessions(s).labels;
        
        conditions = unique(y);
        n_trials = histc(y, conditions);

        for f = 1:n_perm
            fprintf('\tPermutation %g / %g\n', f, n_perm)
            
            % precompute permutations
            ind_pseudo_train = nan(n_conditions, n_conditions, 2*(n_pseudo-1));
            ind_pseudo_test = nan(n_conditions, n_conditions, 2);
            labels_pseudo_train = nan(n_conditions, n_conditions, 2*(n_pseudo-1));
            labels_pseudo_test = nan(n_conditions, n_conditions, 2);
            for c1 = 1:n_conditions
                range_c1 = (c1-1)*(n_pseudo-1)+1:c1*(n_pseudo-1);
                for c2 = 1:n_conditions
                    range_c2 = (c2-1)*(n_pseudo-1)+1:c2*(n_pseudo-1);
                    ind_pseudo_train(c1, c2, 1:2*(n_pseudo - 1)) = [range_c1 range_c2];
                    ind_pseudo_test(c1, c2, :) = [c1 c2];
                    labels_pseudo_train(c1, c2, 1:2*(n_pseudo - 1)) = ...
                        [conditions(c1)*ones(1,n_pseudo-1) conditions(c2)*ones(1,n_pseudo-1)];
                    labels_pseudo_test(c1, c2, :) = conditions([c1 c2]);
                end
            end              
            train_indices = cell(1, n_conditions*(n_pseudo-1));
            test_indices = cell(1, n_conditions);
            for c1 = 1:n_conditions  % separate permutation for each condition
                prm_ = randperm(n_trials(c1));                
                prm = cell(1, n_pseudo);
                splitsize = n_trials(c1) / n_pseudo;
                for i = 1:n_pseudo
                    idxs = floor(round((i-1)*splitsize)):floor(round((i)*splitsize))-1;
                    prm{i} = prm_(idxs + 1);
                end                                
                ind = cellfun(@(x)x+sum(n_trials(1:c1-1)), prm, 'UniformOutput', 0);
                xrange = (c1-1)*(n_pseudo-1)+1:c1*(n_pseudo-1);
                for i = 1:length(xrange)
                    train_indices{xrange(i)} = ind{i}; % train indices: 80% aller trials beider Bedigungen pro run: die ersten 4 trials sind aus Bedingung A, die letzten 4 aus Bedingung B(jeweils 80%)
                end
                test_indices{c1} = ind{end}; % 20 % der trials beider Bedingungen gemischt
            end                                

            % 1. Compute pseudo-trials for training and test
            Xpseudo_train = nan(length(train_indices), n_sensors, n_time);
            Xpseudo_test = nan(length(test_indices), n_sensors, n_time);
            for i = 1:length(train_indices)
                Xpseudo_train(i, :, :) = mean(X(train_indices{i}, :, :), 1);
            end
            for i = 1:length(test_indices)
                Xpseudo_test(i, :, :) = mean(X(test_indices{i}, :, :), 1);
            end


            % 2. Whitening using the Epoch method
            sigma_conditions = reshape(squeeze(labels_pseudo_train(1,:,n_pseudo:end))',1,[]);
            sigma_ = nan(n_conditions, n_sensors, n_sensors);
            for c = 1:n_conditions
                % compute sigma for each time point, then average across time
                tmp_ = nan(n_time, n_sensors, n_sensors);
                for t = 1:n_time
                    tmp_(t, :, :) = cov1para(Xpseudo_train(sigma_conditions==c, :, t));
                end
                sigma_(c, :, :) = mean(tmp_, 1);
            end
            sigma = squeeze(mean(sigma_, 1));  % average across conditions
            sigma_inv = sigma^-0.5;
            for t = 1:n_time
                Xpseudo_train(:, :, t) = squeeze(Xpseudo_train(:, :, t)) * sigma_inv;
                Xpseudo_test(:, :, t) = squeeze(Xpseudo_test(:, :, t)) * sigma_inv;
            end

            for t = 1:n_time
                for c1 = 1:n_conditions-1
                    for c2 = c1+1:n_conditions
                        % 3. Fit the classifier using training data
                        data_train = Xpseudo_train(ind_pseudo_train(c1, c2, :), :, t);
                        y_train = squeeze(labels_pseudo_train(c1, c2, :));
                        model_svm = svmtrain(y_train, data_train, '-c 1 -q 0 -t 0'); % libSVM: -c cost : set the parameter C of C-SVC, epsilon-SVR, and nu-SVR (default 1)\n"; -q : quiet mode (no outputs)\; "-t kernel_type : set type of kernel function (default 2)\n"                       
                        model_weird = weirdtrain(y_train, data_train);
                        model_gnb = gnbtrain(y_train, data_train);

                        % 4. Compute and store classification accuracies
                        data_test = Xpseudo_test(ind_pseudo_test(c1, c2, :), :, t);
                        y_train = squeeze(labels_pseudo_test(c1, c2, :));
                        result.svm(s, f, c1, c2, t) = ...
                            mean(svmpredict(y_train,data_test,model_svm,'-q 0 -t 0')==y_train)-0.5;
                        result.weird(s, f, c1, c2, t) = ...
                            mean(weirdpredict(y_train,data_test,model_weird)==y_train)-0.5;
                        result.gnb(s, f, c1, c2, t) = ...
                            mean(gnbpredict(y_train,data_test,model_gnb)==y_train)-0.5;
                    end
                end
            end
        end
    end
    % average across permutations
    for c = 1:length(clfs)
        result_.(clfs{c}) = nan(n_sessions, n_perm, n_conditions, n_conditions, n_time);
    end 
    result_.svm = squeeze(nanmean(result.svm, 2));
    result_.gnb = squeeze(nanmean(result.gnb, 2));
    result_.weird = squeeze(nanmean(result.weird, 2));
    result = result_; % result = 4D matrix; Beispiel: 2 *9*9*111 => 2 = Anzahl sessions, 9 = Anzahl an Bedingungen, 111 = Anzahl an Zeitpunkten
    result.time = sessions.time;
    result.like = '1';
    result.dislike = '2';
    result.dontcare = '3';
    result.atlas = atlas;
    
    fn_outPath = [outPath 'MEG\sourcespace\all_clusters\Guggenmos_decoding_results\'];
    if ~exist (fn_outPath, 'dir')
        mkdir (fn_outPath)
    end
    save([fn_outPath 'result_decoding_' freq], 'result')
    
%% Now we plot the average classification accuracy time course by collapsing across sessions and conditions:
if 3 == length(conditions) 
    figure; axis tight
    hold on
    plot(sessions.time{1,1}, 100*squeeze(result.svm(1,2,:))+50, 'b-', 'linewidth', 1, 'DisplayName','like vs dislike'); %legend('like vs dislike');
    plot(sessions.time{1,1}, 100*squeeze(result.svm(1,3,:))+50, 'r-', 'linewidth', 0.5, 'DisplayName', 'like vs dontcare');
    plot(sessions.time{1,1}, 100*squeeze(result.svm(2,3,:))+50, 'k-', 'linewidth',0.5, 'DisplayName', 'dislike vs dontcare');
    legend('show')
    plot([-0.5 1], [50 50], 'k-')
    xlim([-0.5 1])
    xlabel('Time [s]');
    ylabel('Classification accuracy svm');
    savefig([outPath 'MEG\sourcespace\all_clusters\Guggenmos_decoding_results\decoding_result_' freq '_svm.fig'])
    close
    figure; axis tight
    hold on
    plot(sessions.time{1,1}, 100*squeeze(result.weird(1,2,:))+50, 'b-', 'linewidth', 1)
    plot(sessions.time{1,1}, 100*squeeze(result.weird(1,3,:))+50, 'r-', 'linewidth', 0.5)
    plot(sessions.time{1,1}, 100*squeeze(result.weird(2,3,:))+50, 'k-', 'linewidth', 0.5)
    plot([-0.5 1], [50 50], 'k-')
    xlim([-0.5 1])
    xlabel('Time [s]')
    ylabel('Classification accuracy weird')
    legend('like vs dislike', 'like vs dontcare', 'dislike vs dontcare')
    savefig([outPath 'MEG\sourcespace\all_clusters\Guggenmos_decoding_results\decoding_result_' freq '_weird.fig'])
    close
    figure; axis tight
    hold on
    plot(sessions.time{1,1}, 100*squeeze(result.gnb(1,2,:))+50, 'b-', 'linewidth', 1)
    plot(sessions.time{1,1}, 100*squeeze(result.gnb(1,3,:))+50, 'r-', 'linewidth', 0.5)
    plot(sessions.time{1,1}, 100*squeeze(result.gnb(2,3,:))+50, 'k-', 'linewidth', 0.5)
    plot([-0.5 1], [50 50], 'k-')
    xlim([-0.5 1])
    xlabel('Time [s]')
    ylabel('Classification accuracy gnb')
    legend('like vs dislike', 'like vs dontcare', 'dislike vs dontcare')
    savefig([outPath 'MEG\sourcespace\all_clusters\Guggenmos_decoding_results\decoding_result_' freq '_gnb.fig'])
    close
else
    figure; axis tight
    hold on
    plot(sessions.time{1,1}, 100*squeeze(nanmean(nanmean(result.svm, 1), 2))+50, 'linewidth', 1) % nanmean = mean ignoring nans
    plot(sessions.time{1,1}, 100*squeeze(nanmean(nanmean(result.weird, 1), 2))+50, 'linewidth', 1)
    plot(sessions.time{1,1}, 100*squeeze(nanmean(nanmean(result.gnb, 1), 2))+50, 'linewidth', 1)
    plot([-0.5 1], [50 50], 'k-')
    xlim([-0.5 1])
    xlabel('Time [s]')
    ylabel('Classification accuracy')
    legend('SVM', 'WeiRD', 'GNB')
    savefig([outPath 'MEG\sourcespace\all_clusters\Guggenmos_decoding_results\decoding_result_' freq 'SVM_WeiRD_BNB.fig'])
end

% Already for one participant and a reduced data set (10 insteada of 92 conditions), these results look 
% like canonical decoding time courses. Note that the period -100ms to 0ms is the baseline phase and stimulus onset is at 0ms. In this example, Support Vector Machine and WeiRD outperform Gaussian Naive Bayes.

 end
 %%
 

function [vs_allRuns] = adi_appendvirtsensMEG(path2data, outPath, freqname, condition, atlas)

% append Runs 
% Channels   --> Configuration of Brainstorm MEG Channels    
if ~exist([outPath 'MEG\sourcespace\noROIs\runs_appended\virtsens\virtsens_all_conditions_allRuns_', freqname, '.mat'], 'file')
    vs_allRuns = [];
    spatialfilter = {};

    for j = 1:size(condition,2)
        fileName = ([condition{j} '*.mat']);
        files = dir(fullfile(path2data, fileName));
        size_files = size(files);

        for i = 1:size_files(1,1)
            load ([path2data files(i).name]);

            [data_bpfreq] = adi_bpfilter(cleanMEG_interp, freqname);
         
            if ~isfield(spatialfilter, (['run' num2str(i)]))
                spatialfilter.(['run' num2str(i)]) = load ([outPath 'MEG\sourcespace\spatialfilter\run' num2str(i) '\spatialfilter_loose_orientation_singletrials_' freqname '.mat']);
                spatialfilter.(['run' num2str(i)]) = cat(1,spatialfilter.(['run' num2str(i)]).spatialfilter_orig{:});   
            end

                
            virtsens = [];
            for k = 1:length(data_bpfreq.trial)
                virtsens.trial{k} = spatialfilter.(['run' num2str(i)])*data_bpfreq.trial{k};
            end
        
            virtsens.time = data_bpfreq.time;
            virtsens.fsample = data_bpfreq.fsample;
            for p = 1:length(virtsens.trial)
                euclid_norm = zeros(length(virtsens.trial{1,p})/3, size(virtsens.trial{1,p}(1,:),2));
                n = 1;
                for k = 1:3:length(virtsens.trial{1,p})
                    euclid_norm(n,:) = sqrt(virtsens.trial{1,p}(k,:).^2+virtsens.trial{1,p}(k+1,:).^2+virtsens.trial{1,p}(k+2,:).^2);
                    n = n+1;
                end  
                vs_euclid_norm.trial{1,p} = euclid_norm;
                clear euclid_norm
            end
            virtsens.trial = vs_euclid_norm.trial;
            
            virtsens.label = cell(length(virtsens.trial{1,1}(:,1)), 1);
            for n = 1:length(virtsens.trial{1,1}(:,1))
                virtsens.label{n,1}=num2str(n);
            end
            
            % sanity check:
            cfg=[];
            avg = ft_timelockanalysis(cfg, virtsens);
            figure;
            plot(virtsens.time{1,1}, mean(avg.avg))
            sanity_path = [outPath 'MEG\sourcespace\all_voxels_without_cerebellum\sanity_check\'];
            if ~exist(sanity_path, 'dir')
                mkdir (sanity_path)
            end
            savefig([sanity_path 'mean_virtsens_' files(i).name '_' freqname '.fig'])
            close
            
         % siehe Skript von Yuval: => all will have a similar noise level
        
            if 1 == strcmp (cfg_virtsens, 'virtsens_ns')
                ns = mean(abs(spatialfilter),2);

                for k = 1:length(virtsens.trial)
                    virtsens_ns.trial{k} = virtsens.trial{1,k}./repmat(ns,1,size(virtsens.trial{1,k},2)); % => all will have a similar noise level
                end

                virtsens_ns.time = virtsens.time;
                virtsens_ns.fsample = virtsens.fsample;

                for k = 1:length(virtsens_ns.trial{1}(:,1))
                    virtsens_ns.label{k} = num2str(k);
                end
                virtsens = virtsens_ns;
                clear virtsens_ns
            end

            files2append.(condition{j}).(['run' num2str(i)]) = virtsens;
            clear virtsens data_bpfreq cleanMEG_interp vs_euclid_norm

        end
    end 
    condition = fields(files2append);
    
    for k = 1:size(fields(files2append),1)
        switch size(fields(files2append.(condition{k})),1)
            case 3 
                vs_allRuns.(condition{k}).trial = [files2append.(condition{k}).(['run' num2str(1)]).trial files2append.(condition{k}).(['run' num2str(2)]).trial files2append.(condition{k}).(['run' num2str(3)]).trial];
                vs_allRuns.(condition{k}).time = [files2append.(condition{k}).(['run' num2str(1)]).time files2append.(condition{k}).(['run' num2str(2)]).time files2append.(condition{k}).(['run' num2str(3)]).time];
                vs_allRuns.(condition{k}).label =  files2append.(condition{k}).(['run' num2str(1)]).label;
            case 2
                vs_allRuns.(condition{k}).trial = [files2append.(condition{k}).(['run' num2str(1)]).trial files2append.(condition{k}).(['run' num2str(2)]).trial];
                vs_allRuns.(condition{k}).time = [files2append.(condition{k}).(['run' num2str(1)]).time files2append.(condition{k}).(['run' num2str(2)]).time];
                vs_allRuns.(condition{k}).label =  files2append.(condition{k}).(['run' num2str(1)]).label;
        end
    end
    
    % sanity check Nr. 2:
   
    cfg = [];
    for m = 1:size(fields(vs_allRuns),1)
        avg = ft_timelockanalysis(cfg, vs_allRuns.(condition{m}));
        figure;
        plot(avg.time, mean(avg.avg))
        savefig([sanity_path 'mean_virtsens_allRuns_' (condition{m}) '_' freqname '.fig'])
        close
    end
    
    pathAppended = [outPath '\MEG\sourcespace\all_voxels_without_cerebellum\runs_appended\virtsens\'];
    if ~exist ( pathAppended, 'dir')
        mkdir (pathAppended)
    end
   
    save ([outPath 'MEG\sourcespace\all_voxels_without_cerebellum\runs_appended\virtsens\' cfg_virtsens '_all_conditions_allRuns_', freqname, '.mat'], 'vs_allRuns', '-v7.3');

else
    vs_allRuns = [];
    %load ([outPath 'MEG\sourcespace\all_voxels_without_cerebellum\runs_appended\virtsens\' cfg_virtsens '_all_conditions_allRuns_', freqname, '.mat'], 'vs_allRuns');
end


% session.data = session.data(:, 92:end,:);
% session.atlas = atlas;
end





function [data_bpfreq_res_sel] = adi_bpfilter(filename, bpname)


switch bpname
    case 'bp1-45Hz'
        bpfreq = 45;
    case 'delta'
        bpfreq = 4;
    case 'theta'
        bpfreq = [4 8];
    case 'alpha'
        bpfreq = [8 13];
    case 'beta'
        bpfreq = [13 25];
    case 'low_gamma'
        bpfreq = [25 45];
    case 'high_gamma'
        bpfreq = [55 90];
end

cfg=[];
cfg.trials        = 'all'; 
cfg.feedback = 'yes';
if 1 == strcmp(bpname, 'delta') || 1 == strcmp(bpname, 'bp1-45Hz')
    cfg.lpfilter      = 'yes';
    cfg.lpfreq        = bpfreq;
else
    cfg.bpfilter      = 'yes'; 
    cfg.bpfreq        = bpfreq;
end

[data_bpfreq] = ft_preprocessing(cfg, filename);
data_bpfreq.ChannelFlag_Bst = filename.ChannelFlag_Bst;
if isfield(filename, 'trialinfo')
    data_bpfreq.trialinfo = filename.trialinfo;
end

for k = 1:length(data_bpfreq.trial)
    data_bpfreq.trial{1,k} = data_bpfreq.trial{1,k}(1:248,:);
end    

for k = 1:length(data_bpfreq.trial)
    data_bpfreq.grad.label(249:end) = [];
    data_bpfreq.grad.chanori(249:end, :) = [];
    data_bpfreq.grad.chanpos(249:end, :) = [];
    data_bpfreq.grad.tra(249:end, :) = [];
    data_bpfreq.label(249:end) = [];
end

for m = 1:length(data_bpfreq.trial)
    if 3053 == length(data_bpfreq.time{1,m}) 
        data_bpfreq.trial{1,m}(:,1) = [];
        data_bpfreq.time{1,m}(:,1) = [];
    end
end
cfg =[];
cfg.resamplefs = 256;
cfg.detrend = 'no';
[data_bpfreq_res] = ft_resampledata(cfg, data_bpfreq);

cfg =[];
cfg.latency = [-0.5 1];
data_bpfreq_res_sel = ft_selectdata(cfg, data_bpfreq_res);

       
end
     
