function [CV] = adi_crossval(like, dislike, method, num_holdout)


switch method
    
    case 'holdout'
        
        %% like
        [count_trl_like_subject, subj_like] = histcounts(categorical(like), categorical(unique(like)));
        for k = 1:length(subj_like)
            subj_like_array(k) = str2num(cell2mat(subj_like(k)));
        end
        [CV.like] = CV_holdout(subj_like_array, num_holdout);
        
        %% dislike
        [count_trl_dislike_subject, subj_dislike] = histcounts(categorical(dislike), categorical(unique(dislike)));
        for k = 1:length(subj_dislike)
            subj_dislike_array(k) = str2num(cell2mat(subj_dislike(k)));
        end
        [CV.dislike] = CV_holdout(subj_dislike_array, num_holdout);
        
         %% w�hle alle m�glichen Kombinationen aus:
         size_datasplit_like = 1:size(CV.like,2);
         size_datasplit_dislike = 1:size(CV.dislike,2);
         [p,q] = meshgrid(size_datasplit_like, size_datasplit_dislike);
         CV.pairs = [p(:) q(:)];
       
        
    case 'holdout_balldesign'   
        
        
        
        
        
        
        
    otherwise
        return
         
end
    
   
end


%% helperfunctions



function [CV] = CV_holdout(subj_array, num_holdout)

    subjects = subj_array;
    [test_rand, indx] = datasample(subjects, num_holdout, 'Replace', false);
    CV(1).test = test_rand;
    test_set = subjects;
    training_set = subjects;
    training_set(find(ismember(training_set, test_rand)))=[]; 
    CV(1).training = training_set;
    test_set(indx)=[];
    i = 1;
    while ~isempty(test_set) && size(test_set,2) > 1 
        i = i + 1;
        [test_rand, indx] = datasample(test_set, num_holdout, 'Replace', false);
        CV(i).test = test_rand;
        training_set = subjects;
        training_set(find(ismember(training_set, test_rand))) = []; 
        CV(i).training = training_set; 
        test_set(indx) = [];
        clear test_rand indx
    end
    
    if 1 == size(test_set,2)
        i = i + 1;
        [test_rand, indx] = datasample(test_set, 1, 'Replace', false);
        CV(i).test = test_rand;
        training_set = subjects;
        training_set(find(ismember(training_set, test_rand))) = [];
        CV(i).training = training_set; 
        test_set(indx) = [];
        clear test_rand indx
    end
        
end