%%


tbl = table(list_regression.Proband, list_regression.favorite_color, list_regression.favorite_pattern, list_regression.Congruency_favorite_Ratings,...
    list_regression.Farbkombi, list_regression.Muster, list_regression.Balldesign, list_regression.rank_fav,list_regression.like_ratio,list_regression.xnmax_condition,...
    list_regression.condition_num, list_regression.perf_comp1,'VariableNames',{'subject', 'favorite_color', 'favorite_pattern', 'Congruency_favorite_Ratings', 'color', 'pattern', 'Balldesign','rank_fav','like_ratio','trialnumber','like_vs_dislike','perf_comp1'});

tbl.subject = categorical(tbl.subject);
tbl.like_vs_dislike = categorical(tbl.like_vs_dislike);
tbl.Balldesign = categorical(tbl.Balldesign);
tbl.pattern = categorical(tbl.pattern);
tbl.color = categorical(tbl.color);
tbl.Congruency_favorite_Ratings = categorical(tbl.Congruency_favorite_Ratings);
tbl.favorite_color = categorical(tbl.favorite_color);
tbl.favorite_pattern = categorical(tbl.favorite_pattern);

% rank_fav hat keinen Einfluss
model = 'perf_comp1~trialnumber+Farbkombi+Muster+Farbkombi:Muster'; 
model = 'perf_comp1~trialnumber+like_ratio'; 
model = 'perf_comp1~trialnumber+like_vs_dislike+like_vs_dislike:Congruency_favorite_Ratings+'; 
model = 'perf_comp1~trialnumber+like_ratio+like_ratio:Congruency_favorite_Ratings+color+pattern+color:pattern';
model = 'perf_comp1~trialnumber+like_ratio+like_ratio:Congruency_favorite_Ratings+color+pattern+color:pattern+color:favorite_color+pattern:favorite_pattern';

lm = fitlm(tbl,model, 'RobustOpts','on')

model = 'perf_comp1~color+pattern+like_ratio+trialnumber+favorite_pattern+favorite_color+color:pattern+favorite_color:color'; % always with constant
lm = fitlm(tbl,model, 'RobustOpts','on')


%% regression favorites:



tbl_fav = table(list_regression.subject, list_regression.congruency_new, list_regression.color, list_regression.pattern, list_regression.no_Trial, ...
    list_regression.like_ratio, list_regression.like_dichotom, list_regression.fav_pattern, list_regression.fav_color, list_regression.fav_rank,  ...
    list_regression.perf_comp2, 'VariableNames',{'subject', 'congruency_fav_ratings', 'color', 'pattern', 'no_Trial',...
    'like_ratio', 'like_vs_dislike', 'fav_pattern', 'fav_color','rank_fav','perf_comp2'});


tbl_fav.subject = categorical(tbl_fav.subject);
tbl_fav.like_vs_dislike = categorical(tbl_fav.like_vs_dislike);
tbl_fav.pattern = categorical(tbl_fav.pattern);
tbl_fav.color = categorical(tbl_fav.color);
tbl_fav.congruency_fav_ratings = categorical(tbl_fav.congruency_fav_ratings);
tbl_fav.fav_color = categorical(tbl_fav.fav_color);
tbl_fav.fav_pattern = categorical(tbl_fav.fav_pattern);
% Stefan: keine 4er-Interaktion, sondern 3er-Interaktion nehmen:
% model = ['perf_comp1~no_Trial+Farbkombi+Muster+Farbkombi:Muster+rank_fav+rank_fav:congruency_fav_ratings+Farbkombi:Muster:rank_fav+Farbkombi:Muster:rank_fav:congruency_fav_ratings ']; 
model = ['perf_comp2~no_Trial+rank_fav+color+pattern+color:pattern+color:pattern:rank_fav ']; 
model = ['perf_comp2~no_Trial+like_vs_dislike+color+pattern+color:pattern+rank_fav:like_vs_dislike+color:pattern:rank_fav '];


model = ['perf_comp2~no_Trial+color+pattern+color:pattern+rank_fav+rank_fav:congruency_fav_ratings+color:pattern:rank_fav ']; % best

model = ['perf_comp2~no_Trial+color+pattern+color:pattern+rank_fav+rank_fav:congruency_fav_ratings+color:pattern:rank_fav ']; % best
% model = ['perf_comp2~no_Trial+color+pattern+color:pattern+rank_fav+color:pattern:rank_fav ']; 
model = ['perf_comp2~no_Trial+like_vs_dislike']
lm_tst = fitlm(tbl_fav,model, 'weights', [0.5*ones(size(tbl_fav,1)/2,1); ones(size(tbl_fav,1)/2,1)  ])
lm = fitlm(tbl_fav,model, 'Weights', [ones(size(tbl_fav,1)/2,1); ones(size(tbl_fav,1)/2,1)  ])



model = ['perf_comp2~no_Trial+rank_fav+color+pattern+color:pattern+rank_fav:congruency_fav_ratings+color:pattern:rank_fav ']; 
lm = fitlm(tbl_fav,model, 'RobustOpts','off') % like = dichotom entspricht alter congruency

lm = fitlm(tbl_fav,model, 'RobustOpts','on')


boxplot(tbl_fav.perf_comp2, tbl_fav.like_vs_dislike)
find(list_regression.like_dichotom == 2)


boxplot(tbl_fav.perf_comp2, tbl_fav.rank_fav)










scatter(tbl.trialnumber, tbl.perf_comp1)


%%

X=[list_regression.rank_fav list_regression.like_ratio list_regression.xnmax_condition];
mdl = fitlm(X,y)
y = list_regression.perf_2xNmax;
% stepwise sollte nicht verwendet werden, intercept sollte immer dabei sein
mdl1 = stepwiselm(tbl,model)
mdl_comp1 = stepwiselm(tbl,'interactions','ResponseVar','perf_comp1'); % best
mdl_comp1 = stepwiselm(tbl,'interactions','ResponseVar','perf_comp1', 'Intercept', 0); % best
% dummy Variablen:
temp_congr_favorite_pattern = dummyvar(categorical(tbl.favorite_pattern));
pattern_congruent = temp_congr_favorite_pattern(:,1);
pattern_incongruent = temp_congr_favorite_pattern(:,2);
tbl.pattern_congruent = pattern_congruent;
tbl.pattern_incongruent = pattern_incongruent;

temp_congr_favorite_color = dummyvar(categorical(tbl.favorite_color));
color_congruent = temp_congr_favorite_color(:,1);
color_incongruent = temp_congr_favorite_color(:,2);
tbl.color_congruent = color_congruent;
tbl.color_incongruent = color_incongruent;



mdl.ObservationInfo(84,:)
temp_Year = dummyvar(categorical(Model_Year));
Model_Year_70 = temp_Year(:,1);
Model_Year_76 = temp_Year(:,2);
Model_Year_82 = temp_Year(:,3);
tbl = table(Model_Year_70,Model_Year_76,Model_Year_82,MPG);
mdl = fitlm(tbl,'MPG ~ Model_Year_70 + Model_Year_76 + Model_Year_82 - 1')


lm = fitlm(tbl,'perf_2xNmax~like_ratio+trialnumber+Farbkombi+Muster+Congruency_favorite_Ratings+favorite_pattern+favorite_color', 'RobustOpts','on')

lm = fitlm(tbl,'perf_comp1~subject+like_ratio+trialnumber+favorite_pattern - 1', 'RobustOpts','on') % best: ohne intercept: r²=.44
lm = fitlm(tbl,'perf_comp1~subject+Balldesign+like_vs_dislike+trialnumber+favorite_pattern - 1', 'RobustOpts','on') 
lm = fitlm(tbl,'perf_2xNmax~like_ratio+trialnumber+pattern_congruent', 'RobustOpts','on') %r²=0.05
lm = fitlm(tbl,'perf_2xNmax~like_ratio+trialnumber+favorite_pattern+favorite_color - 1', 'RobustOpts','on') % r²=.39
lm = fitlm(tbl,'perf_2xNmax~like_ratio+trialnumber+pattern_congruent+color_congruent - 1', 'RobustOpts','on') % %r²=.44


model = 'perf_comp1~subject+like_ratio+trialnumber+favorite_pattern'; %best
model = 'perf_comp1~favorite_color+favorite_pattern+trialnumber'; %best

lm = fitlm(tbl,'perf_comp1~subject+like_ratio+trialnumber+favorite_pattern - 1', 'RobustOpts','on') 


Acceleration:Weight + Weight^2
anova(mdl,'summary')
plotResiduals(mdl2_comp1)
outlier = mdl2_comp1.Residuals.Raw < -0.04;
mdl2_comp1 = fitlm(patients,modelspec,...
    'Exclude',116);
mdl_2_comp1 = stepwiselm(tbl,'interactions','ResponseVar','perf_comp1', 'Exclude', 116) % best

figure
plotInteraction(mdl1,'favorite_pattern','rank_fav')

plotInteraction(mdl,var1,var2)

mdl1 = step(mdl,'NSteps',10)

mdl_glm_comp1 = stepwiseglm(tbl)


find(outlier)