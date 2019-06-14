



for i=1:length(data_like.balldesign)
    switch char(data_like.balldesign{i})
        
        case balldesigns_names{1}
            data_like.balldesign_num(i) = 1;
        case balldesigns_names{2}
            data_like.balldesign_num(i) = 2;
        case balldesigns_names{3}
            data_like.balldesign_num(i) = 3;
        case balldesigns_names{4}
            data_like.balldesign_num(i) = 4;
    end
   
end

counts=histcounts(data_like.balldesign_num);


data_like_first.trial=[data_like.trial(1:162) data_like.trial(326:(326+161)) data_like.trial(counts(1)+counts(2)+1:counts(1)+counts(2)+floor(counts(3)/2))]
data_like_first.time=[data_like.time(1:162) data_like.time(326:(326+161)) data_like.time(counts(1)+counts(2)+1:counts(1)+counts(2)+floor(counts(3)/2))]
data_like_first.balldesign=[data_like.balldesign(1:162) data_like.balldesign(326:(326+161)) data_like.balldesign(counts(1)+counts(2)+1:counts(1)+counts(2)+floor(counts(3)/2))]
data_like_first.label=data_like.label;
data_like_first.response_label=[data_like.response_label(1:162) data_like.response_label(326:(326+161)) data_like.response_label(counts(1)+counts(2)+1:counts(1)+counts(2)+floor(counts(3)/2))]
data_like_first.balldesign_num=[data_like.balldesign_num(1:162) data_like.balldesign_num(326:(326+161)) data_like.balldesign_num(counts(1)+counts(2)+1:counts(1)+counts(2)+floor(counts(3)/2))]

data_like_second.trial=[data_like.trial(163:325) data_like.trial(325+floor(325/2)+1:325+floor(325/2)+162+1) data_like.trial(counts(1)+counts(2)+floor(counts(3)/2)+1:size(data_like.trial,2))]
data_like_second.time=[data_like.trial(163:325) data_like.time(325+floor(325/2)+1:325+floor(325/2)+162+1) data_like.time(counts(1)+counts(2)+floor(counts(3)/2)+1:size(data_like.time,2))]
data_like_second.balldesign=[data_like.trial(163:325) data_like.balldesign(325+floor(325/2)+1:325+floor(325/2)+162+1) data_like.balldesign(counts(1)+counts(2)+floor(counts(3)/2)+1:size(data_like.balldesign,2))]
data_like_second.response_label=[data_like.response_label(163:325) data_like.response_label(325+floor(325/2)+1:325+floor(325/2)+162+1) data_like.response_label(counts(1)+counts(2)+floor(counts(3)/2)+1:size(data_like.response_label,2))]
data_like_second.balldesign_num=[data_like.balldesign_num(163:325) data_like.balldesign_num(325+floor(325/2)+1:325+floor(325/2)+162+1) data_like.balldesign_num(counts(1)+counts(2)+floor(counts(3)/2)+1:size(data_like.balldesign_num,2))]
data_like_second.label=data_like.label;


%----------------------------

for i=1:length(data_dislike.balldesign)
    switch char(data_dislike.balldesign{i})
        
        case balldesigns_names{1}
            data_dislike.balldesign_num(i) = 1;
        case balldesigns_names{2}
            data_dislike.balldesign_num(i) = 2;
        case balldesigns_names{3}
            data_dislike.balldesign_num(i) = 3;
        case balldesigns_names{4}
            data_dislike.balldesign_num(i) = 4;
    end
   
end

counts=histcounts(data_dislike.balldesign_num);


data_dislike_first.trial=[data_dislike.trial(1:162) data_dislike.trial(326:(326+161)) data_dislike.trial(counts(1)+counts(2)+1:counts(1)+counts(2)+floor(counts(3)/2))]
data_dislike_first.time=[data_dislike.time(1:162) data_dislike.time(326:(326+161)) data_dislike.time(counts(1)+counts(2)+1:counts(1)+counts(2)+floor(counts(3)/2))]
data_dislike_first.balldesign=[data_dislike.balldesign(1:162) data_dislike.balldesign(326:(326+161)) data_dislike.balldesign(counts(1)+counts(2)+1:counts(1)+counts(2)+floor(counts(3)/2))]
data_dislike_first.label=data_dislike.label;
data_dislike_first.response_label=[data_dislike.response_label(1:162) data_dislike.response_label(326:(326+161)) data_dislike.response_label(counts(1)+counts(2)+1:counts(1)+counts(2)+floor(counts(3)/2))]
data_dislike_first.balldesign_num=[data_dislike.balldesign_num(1:162) data_dislike.balldesign_num(326:(326+161)) data_dislike.balldesign_num(counts(1)+counts(2)+1:counts(1)+counts(2)+floor(counts(3)/2))]

data_dislike_second.trial=[data_dislike.trial(163:325) data_dislike.trial(325+floor(325/2)+1:325+floor(325/2)+162+1) data_dislike.trial(counts(1)+counts(2)+floor(counts(3)/2)+1:size(data_dislike.trial,2))]
data_dislike_second.time=[data_dislike.trial(163:325) data_dislike.time(325+floor(325/2)+1:325+floor(325/2)+162+1) data_dislike.time(counts(1)+counts(2)+floor(counts(3)/2)+1:size(data_dislike.time,2))]
data_dislike_second.balldesign=[data_dislike.trial(163:325) data_dislike.balldesign(325+floor(325/2)+1:325+floor(325/2)+162+1) data_dislike.balldesign(counts(1)+counts(2)+floor(counts(3)/2)+1:size(data_dislike.balldesign,2))]
data_dislike_second.response_label=[data_dislike.response_label(163:325) data_dislike.response_label(325+floor(325/2)+1:325+floor(325/2)+162+1) data_dislike.response_label(counts(1)+counts(2)+floor(counts(3)/2)+1:size(data_dislike.response_label,2))]
data_dislike_second.balldesign_num=[data_dislike.balldesign_num(163:325) data_dislike.balldesign_num(325+floor(325/2)+1:325+floor(325/2)+162+1) data_dislike.balldesign_num(counts(1)+counts(2)+floor(counts(3)/2)+1:size(data_dislike.balldesign_num,2))]
data_dislike_second.label=data_dislike.label;







[indx_min_counts ] = find(counts==min(counts));
trls_min = min(counts);

like=[data_cond.trial(randperm(counts(1), trls_min)) data_cond.trial(counts(1)+randperm(counts(2), trls_min)) data_cond.trial(counts(1)+counts(2)+randperm(counts(3), trls_min))];
balldesign_like=[cond.(condition).balldesign(randperm(counts(1), trls_min)) cond.(condition).balldesign(counts(1)+randperm(counts(2), trls_min)) cond.(condition).balldesign(counts(1)+counts(2)+randperm(counts(3), trls_min)) ];
time_red = [data_cond.time(1:trls_min) data_cond.time(counts(1)+1:counts(1)+trls_min) data_cond.time(counts(1)+counts(2)+1:counts(1)+counts(2)+trls_min)];
response_like=[cond.(condition).response_label(1:trls_min) cond.(condition).response_label(counts(1)+1:counts(1)+trls_min) cond.(condition).response_label(counts(1)+counts(2)+1:counts(1)+counts(2)+trls_min)];

data_like.trial=like;
data_like.balldesign=balldesign_like;
data_like.time = time_red;
data_like.response_label = response_like;
data_like.label=label';

%% interaction effect

data_like1=data_like_first
for k=1:length(data_like_first.trial)
    data_like1.individual(k,:,:) = data_like_first.trial{k};

end


for k=1:length(data_like2.trial)
    data_like2.individual(k,:,:) = data_like2.trial{k};

end

%Copy GA11 to GAdiff11_12 and perform the assignment GAdiff11_12.individual=GA11.individual-GA12.individual.
data_like2.individual=data_like2.individual(1:486,:,:);
GAdiff11_12 = data_like1;
GAdiff11_12.individual=data_like1.individual-data_like2.individual;


%%

for k=1:length(data_dislike1.trial)
    data_dislike1.individual(k,:,:) = data_dislike1.trial{k};

end


for k=1:length(data_dislike2.trial)
    data_dislike2.individual(k,:,:) = data_dislike2.trial{k};

end

%Copy GA11 to GAdiff11_12 and perform the assignment GAdiff11_12.individual=GA11.individual-GA12.individual.
data_dislike2.individual=data_dislike2.individual(1:486,:,:);
GAdiff21_22 = data_dislike1;
GAdiff21_22.individual=data_dislike1.individual-data_dislike2.individual;

figure
plot(GAdiff21_22.time{1,1}, squeeze(mean(rms(GAdiff21_22.individual))))
hold on
plot(GAdiff11_12.time{1,1}, squeeze(mean(rms(GAdiff11_12.individual))))
legend({'diff_like'; 'diff_dislike'})

figure
plot(data_like1.time{1,1}, mean(squeeze(rms(data_like1.individual))))
hold on
plot(data_like1.time{1,1}, mean(squeeze(rms(data_like2.individual))))
legend({'data_like1'; 'data_like2'})

figure
plot(data_dislike1.time{1,1}, mean(squeeze(rms(data_dislike1.individual))))
hold on
plot(data_dislike1.time{1,1}, mean(squeeze(rms(data_dislike2.individual))))
legend({'data_dislike1'; 'data_dislike2'})