% nicht verwendet, kann gelöscht werden
function [] = sensor_single_ball()


switch balldesign
    
    case 'gbf'
        delete_subjects = {'nl_adi_04'; 'nl_adi_05'; 'nl_adi_10'; 'nl_adi_19';  'nl_adi_28'};
        for i=1:length(subject_list)
            ind(i) = sum(strcmp (subject_list(i).name, delete_subjects));
        end
        subject_list(find(ind),:) = [];
        clear ind
        delete_runs={'nl_adi_17', 1; 'nl_adi_21', [1 2]; 'nl_adi_25', 1; 'nl_adi_29', [2 3]};

    case 'rwv'
        delete_subjects = {'nl_adi_05'};
        for i=1:length(subject_list)
            ind(i) = sum(strcmp (subject_list(i).name, delete_subjects));
        end
        subject_list(find(ind),:) = [];
        clear ind
        delete_runs={'nl_adi_04', 3; 'nl_adi_10', 1; 'nl_adi_13', 1; 'nl_adi_14' 1; 'nl_adi_15', [1 2]; 'nl_adi_17', [1]; 'nl_adi_19', [1]; 'nl_adi_22', 2; 'nl_adi_24', [1 2]; 'nl_adi_25', [1]; 'nl_adi_28', [1]; 'nl_adi_30', [1 3];} ; 


    case 'rws'
        delete_subjects = {'nl_adi_20'; 'nl_adi_29'; 'nl_adi_34'};
        for i=1:length(subject_list)
            ind(i) = sum(strcmp (subject_list(i).name, delete_subjects));
        end
        subject_list(find(ind),:) = [];
        clear ind

        delete_runs={'nl_adi_04', 3; 'nl_adi_05', 1; 'nl_adi_10', 1; 'nl_adi_17' 1; 'nl_adi_22', [2 3]; 'nl_adi_23', [1 2]; 'nl_adi_28', [1]; 'nl_adi_33', 1} ; 

    case 'ggv'
        delete_subjects = {'nl_adi_07'; 'nl_adi_22'; 'nl_adi_34'};
        for i=1:length(subject_list)
            ind(i) = sum(strcmp (subject_list(i).name, delete_subjects));
        end
        subject_list(find(ind),:) = [];
        clear ind

        delete_runs={'nl_adi_04', 3; 'nl_adi_13', 1; 'nl_adi_14', 1; 'nl_adi_15' [1 2]; 'nl_adi_17', [1]; 'nl_adi_20', [1 2]; 'nl_adi_25', [1]} ; 

    case 'gbs'
        delete_subjects = {'nl_adi_31'};
        for i=1:length(subject_list)
            ind(i) = sum(strcmp (subject_list(i).name, delete_subjects));
        end
        subject_list(find(ind),:) = [];
        clear ind

        delete_runs={'nl_adi_05', 1; 'nl_adi_10', 1; 'nl_adi_15', [1 2]; 'nl_adi_21', 1} ; 

    case 'gbv'
        delete_subjects = {'nl_adi_05'; 'nl_adi_10'; 'nl_adi_15'; 'nl_adi_22'};
        for i=1:length(subject_list)
            ind(i) = sum(strcmp (subject_list(i).name, delete_subjects));
        end
        subject_list(find(ind),:) = [];
        clear ind
        delete_runs={'nl_adi_07', 1; 'nl_adi_09', 3; 'nl_adi_13', 1; 'nl_adi_17', 1; 'nl_adi_25', 1; 'nl_adi_33', 1; 'nl_adi_34', 1} ; 


        
    
    case 'rwf'
        subject_list([2 16 17],:) = [];
        delete_subjects = {};
        for i=1:length(subject_list)
            ind(i) = sum(strcmp (subject_list(i).name, delete_subjects));
        end
        subject_list(find(ind),:) = [];
        clear ind
        delete_runs={}; 

    case 'ggf'
        delete_subjects = {'nl_adi_04'; 'nl_adi_07'; 'nl_adi_12'; 'nl_adi_18'; 'nl_adi_20'};
        for i=1:length(subject_list)
            ind(i) = sum(strcmp (subject_list(i).name, delete_subjects));
        end
        subject_list(find(ind),:) = [];
        clear ind

        delete_runs={'nl_adi_17', 1; 'nl_adi_19', 1;  'nl_adi_21', [2 3]; 'nl_adi_23', 3; 'nl_adi_28' [1 2]; 'nl_adi_29', 1; 'nl_adi_34', [2 3]}; % adi28: run2 sollte auch gelöscht werden, existiert aber nicht in den Daten

    case 'ggs';
        delete_subjects = {'nl_adi_20'; 'nl_adi_29'; 'nl_adi_34'};
        for i=1:length(subject_list)
            ind(i) = sum(strcmp (subject_list(i).name, delete_subjects));
        end
        subject_list(find(ind),:) = [];
        clear ind

        delete_runs={'nl_adi_04', 3; 'nl_adi_05', [2];  'nl_adi_10', [1]; 'nl_adi_15', 2; 'nl_adi_17' 1; 'nl_adi_22', [2 3]; 'nl_adi_23', [1 2]; 'nl_adi_28', [1]; 'nl_adi_33', 1} ; 
end






end
