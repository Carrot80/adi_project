%% additional interpolation
clear
subjectpath = ['W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\nl_adi_23\MEG_analysis\noisereduced\1_95Hz\02_interpolated\'];
file_list = dir([subjectpath '*.mat']);

for kk = 1%:length(file_list)
   load([file_list(kk).folder filesep file_list(kk).name]) 
   [cleanMEG_interp] = adi_additionalMEG_interpolation(cleanMEG_interp, []);
   save([file_list(kk).folder filesep file_list(kk).name], 'cleanMEG_interp')
end



%% 24.7. crossvalidation leave out one exemplar per subject with feature reduction:

clear
path2subjects = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';

subjects_dir = dir(path2subjects);
subjects_dir(1:2,:) = [];
balldesign = {'gbf'; 'gbs'; 'gbv'; 'ggf'; 'rwv'; 'rws'; 'rwf'; 'ggs'; 'ggv'};
% adi_leave_out_exemplar_singlsubj(path2subjects, subjects_dir, balldesign)
adi_leave_out_exemplar_singlsubj_ohne_pca(path2subjects, subjects_dir, balldesign)



%% 
figure
plot(1:385, perf.lda.accuracy(4,:))

figure
plot(session.time{1,1}, perf.lda.accuracy(1:3,:))

figure
plot(session.time{1,1}, perf.lda.accuracy(1:3,:))
legend('gbf', 'gbs', 'gbv')
figure
plot(session.time{1,1}, perf.lda.accuracy(4:6,:))
legend('ggf', 'ggs', 'ggv')
figure
plot(session.time{1,1}, perf.lda.accuracy(7:9,:))
legend('rwf', 'rws', 'rwv')

figure
plot(session.time{1,1}, perf.lda.accuracy([1, 4, 7],:))

figure
plot(session.time{1,1}, perf.lda.accuracy([2, 5, 8],:))

figure
plot(session.time{1,1}, perf.lda.accuracy([3, 6, 9],:))

balldesign = {'gbf', 'gbs', 'gbv', 'ggf', 'ggs', 'ggv', 'rwf', 'rws', 'rwv'};
design = {'f', 's', 'v', 'f', 's', 'v', 'f', 's', 'v'};
color = {'gb', 'gb', 'gb', 'gg', 'gg', 'gg', 'rw', 'rw', 'rw'};
lda = perf.lda.accuracy(:, 140:304);

T_acc = table(balldesign',design', color', lda, 'VariableNames',{'balldesign_single','design', 'color', 'accuracy'});
figure
boxplot(lda', balldesign)


%% grand avg all trials
path2subj = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';
subject_list = dir (path2subj);
subject_list([1 2],:) = [];
grandavg = []; 
[grandavg] = adi_grandavg_sensorspace(subject_list, grandavg, []);
save('W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\group_analysis\sensor_space\MEG\grandavg_all_trls\grandavg_all_trls.mat', 'grandavg')

%% grand avg per subject:
path2subj = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';
subject_list = dir (path2subj);
subject_list([1 2],:) = [];
grandavg = []; 
[grandavg] = adi_grandavg_subject(subject_list, grandavg, []);
% save('W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\group_analysis\sensor_space\MEG\grandavg_all_trls\grandavg_all_trls.mat', 'grandavg', '-v7.3')


%% grand avg
path2subj = 'W:\neurochirurgie\science\Kirsten\adidas\fieldtrip_Auswertung\Studie_1_visuell\single_subjects\';
subject_list = dir (path2subj);
subject_list([1 2],:) = [];
folderpath = 'MEG_analysis\noisereduced\1_95Hz\grandavg\';
adi_grandavg_group(subject_list, folderpath);




