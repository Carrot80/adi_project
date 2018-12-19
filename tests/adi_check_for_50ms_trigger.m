function [] = adi_check_for_trigger(path2cleanfile, trigger, subject)
disp(subject)
fileList = dir(fullfile([path2cleanfile, '*.mat']));

for i = 1:length(fileList)
   
    load ([path2cleanfile fileList(i).name])
    disp(fileList(i).name)
    false_trigger = []; %zeros(1, length(cleanMEG.trialinfo.triggerlabel));
    for j=1:length(cleanMEG.trialinfo.triggerlabel)
        ind = find(~ismember(cleanMEG.trialinfo.triggerlabel(j), trigger.triggercodes));
        if ~isempty(ind)
            false_trigger(end+1) = j;
        end
    end
    if any(false_trigger)
        for p = 1:length(false_trigger)
            Msg = ['subject: ' subject ' file: ' fileList(i).name 'trial nr: '  num2str(false_trigger(p)) 'triggercode no ' num2str(cleanMEG.trialinfo.triggerlabel(j)) ];
            fid = fopen(fullfile('W:\neurochirurgie\science\Kirsten\adidas\', 'wrong_triggercodes_new_Data.txt'), 'a');
%             if fid == -1
%               error('Cannot open log file.');
%             end
            fprintf(fid, '%s: %s\n', datestr(now, 0), Msg);
            fclose(fid);
        end
    end
    clearvars cleanMEG
    
end