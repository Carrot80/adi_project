function [] = adi_anova (path, filename)
% see tutorial at bottom
load([path filename])

if 1 == strcmp(filename, 'list_regression_favorites.mat')
    
    like_dichotom_table = list_regression.like_dichotom;
    like_ratio_table = list_regression.like_ratio;
    subject_table = list_regression.Proband;
%     num_trls_table = list_regression.xnmax_condition;
    pattern_table = list_regression.Muster;
    balldesign_table = list_regression.Balldesign;
    farbe_table = list_regression.Farbkombi;
    fav_table = list_regression.rank_fav;
    y_table = list_regression.perf_fav_60_120ms;
    N_table = length(y_table);
  
    figure
    boxplot(y_table, subject_table)  



    
    %% since we have more than 1 measurement per subject: repeated measurement; we should add the subject number as another factor to our n-way anova
%and set it as random factor.
    
    t=table(list_regression_perf.Proband, fav, perf,'VariableNames',{'subject', 'fav','perf'});
m = table([1 2],'VariableNames',{'Measurements'});
m = [1 2]';
rm = fitrm(t,'fav-perf ~ subject',m, 'WithinModel', 'orthogonalcontrasts')
rm = fitrm(t,modelspec)

[p,table,stats] = anovan(list_regression_perf.perf_2xNmax,{list_regression_perf.Balldesign,list_regression_perf.Proband},...
 'random',2,'varnames', {'Balldesign','Subject'});
    
%%   
Perf_comp1 = rand(28,9*3*3*2); % 28 subjects, 9 balldesigns*3pattern*3color*2conditions(likevsdislike)
varNames = cell(9*3*3*2,1);
for i = 1 : 9*3*3*2
 v = strcat('V',num2str(i));
 varNames{i,1} = v;
end

% Create a table storing the respones:
tPerf_comp1 = array2table(Perf_comp1, 'VariableNames',varNames);    
%Create a table reflecting the within subject factors
balldesign = cell(9*3*3*2,1); % 9 balldesigns
pattern = cell(9*3*3*2,1); % 3 patterns
color = cell(9*3*3*2,1); % 3 color combinations    
like =  cell(9*3*3*2,1);

% Assiging the values to the parameters based on the data sorting
c1 = cell(1,1); c1{1} = 'gbf'; c1 = repmat(c1,length(balldesign)/9,1); balldesign(1:length(balldesign)/9*1,1) = c1;
c1 = cell(1,1); c1{1} = 'gbs'; c1 = repmat(c1,length(balldesign)/9,1); balldesign(length(balldesign)/9*1+1:length(balldesign)/9*2,1) = c1;
c1 = cell(1,1); c1{1} = 'gbv'; c1 = repmat(c1,length(balldesign)/9,1); balldesign(length(balldesign)/9*2+1:length(balldesign)/9*3,1) = c1;
c1 = cell(1,1); c1{1} = 'ggf'; c1 = repmat(c1,length(balldesign)/9,1); balldesign(length(balldesign)/9*3+1:length(balldesign)/9*4,1) = c1;
c1 = cell(1,1); c1{1} = 'ggs'; c1 = repmat(c1,length(balldesign)/9,1); balldesign(length(balldesign)/9*4+1:length(balldesign)/9*5,1) = c1;
c1 = cell(1,1); c1{1} = 'ggv'; c1 = repmat(c1,length(balldesign)/9,1); balldesign(length(balldesign)/9*5+1:length(balldesign)/9*6,1) = c1;
c1 = cell(1,1); c1{1} = 'rwf'; c1 = repmat(c1,length(balldesign)/9,1); balldesign(length(balldesign)/9*6+1:length(balldesign)/9*7,1) = c1;
c1 = cell(1,1); c1{1} = 'rws'; c1 = repmat(c1,length(balldesign)/9,1); balldesign(length(balldesign)/9*7+1:length(balldesign)/9*8,1) = c1;
c1 = cell(1,1); c1{1} = 'rws'; c1 = repmat(c1,length(balldesign)/9,1); balldesign(length(balldesign)/9*8+1:length(balldesign)/9*9,1) = c1;

c1 = cell(1,1); c1{1} = 'gelb_blau'; c1 = repmat(c1,length(color)/3,1); color(1:length(color)/3,1) = c1;
c1 = cell(1,1); c1{1} = 'grau_gruen'; c1 = repmat(c1,length(color)/3,1); color(length(color)/3*1+1:length(color)/3*2,1) = c1;
c1 = cell(1,1); c1{1} = 'rot_weiss'; c1 = repmat(c1,length(color)/3,1); color(length(color)/3*2+1:length(color)/3*3,1) = c1;

c1 = cell(1,1); c1{1} = 'Fussball'; c1 = repmat(c1,length(pattern)/9,1); 
c2 = cell(1,1); c2{1} = 'Space'; c2 = repmat(c2,length(pattern)/9,1); 
c3 = cell(1,1); c3{1} = 'Volley'; c3 = repmat(c3,length(pattern)/9,1); 

pattern = [c1; c2; c3; c1; c2; c3; c1; c2; c3];

% like:


    
else
    
like_dichotom_table = list_regression.condition_num;
like_ratio_table = list_regression.like_ratio;
subject_table = list_regression.Proband;
num_trls_table = list_regression.xnmax_condition;
pattern_table = list_regression.Muster;
balldesign_table = list_regression.Balldesign;
farbe_table = list_regression.Farbkombi;
fav_table = list_regression.rank_fav;
y_table = list_regression.perf_60_100ms2xNmax;
N_table = length(y_table);
condition_num = list_regression.condition_num;

figure
boxplot(y_table, subject_table)

%% lösche Einträge mit NaNs in y_table
ind_isnan = find(isnan(y_table));

like_dichotom = like_dichotom_table;
like_dichotom(ind_isnan) = [];
like_ratio = like_ratio_table;
like_ratio(ind_isnan) = [];
subject = subject_table;
subject(ind_isnan) = [];
num_trls = num_trls_table;
num_trls(ind_isnan) = [];
pattern = pattern_table;
pattern(ind_isnan) = [];
balldesign = balldesign_table;
balldesign(ind_isnan) = [];
farbe = farbe_table;
farbe(ind_isnan) = [];
fav = fav_table;
fav(ind_isnan) = [];
y = y_table;
y(ind_isnan) = [];


%% like dislike

ind_dontcare = find(like_dichotom == 3);

like_dichotom_without_dontcare = like_dichotom;
like_dichotom_without_dontcare(ind_dontcare) = [];

pattern_like_dislike = pattern;
pattern_like_dislike(ind_dontcare) = [];

farbe_like_dislike = farbe;
farbe_like_dislike(ind_dontcare) = [];

fav_like_dislike = fav;
fav_like_dislike(ind_dontcare) = [];

y_like_dislike = y;
y_like_dislike(ind_dontcare) = [];

balldesign_like_dislike = balldesign;
balldesign_like_dislike(ind_dontcare) = [];

[p,tbl,stats] = anovan(y_like_dislike, {like_dichotom_without_dontcare }  );
[p,tbl,stats] = anovan(y_like_dislike, {like_dichotom_without_dontcare, balldesign_like_dislike }  );
[p,tbl,stats] = anovan(y_like_dislike, {like_dichotom_without_dontcare, fav_like_dislike }  );
[p,tbl,stats] = anovan(y_like_dislike, {balldesign_like_dislike }  );

% interactions
[p,tbl,stats] = anovan(y_like_dislike,{like_dichotom_without_dontcare farbe_like_dislike pattern_like_dislike},'model','interaction','varnames',{'like_dichotom_without_dontcare','farbe_like_dislike','pattern_like_dislike'})

results = multcompare(stats)

%%
ind_like = find(like_dichotom_without_dontcare == 1);
y_dislike = y_like_dislike;
y_dislike(ind_like)=[];
mean(y_dislike)
std(y_dislike)

ind_dislike = find(like_dichotom_without_dontcare == 2);
y_like = y_like_dislike;
y_like(ind_dislike)=[]; mean(y_like)
std(y_like)

farbe_like = farbe_like_dislike;
farbe_like(ind_dislike)=[];

farbe_dislike = farbe_like_dislike;
farbe_dislike(ind_like)=[];

bar(farbe_dislike)

end

%% alt

figure
scatter(mean_class,like_dislike_ratio)
xlabel('accuracy')
scatter3(mean_class,like_dislike_ratio, fav, 'filled')

% Compute the regression coefficients for a linear model with an interaction term.

X = [ones(size(x1)) x1 x2 x1.*x2];



figure; boxplot(y, fav )
ylabel('mean classification between 30 and 150 ms')
xlabel('rang of favorites')

figure; plot(y, like_ratio, '.b')
axis equal tight
xlabel('mean classification between 30 and 150 ms')
ylabel('like_ratio')

figure; plot(y, num_trls, '.b')
ylabel('number of trials')
xlabel('mean classification between 30 and 150 ms')

% fit to f(x) = b+m*x (m = slope, b = y-intercept):
X = [ones(N,1) like_ratio(:) ];
y = listregression.perf_2n;
subject(find(isnan(y)),:) = [];
X(find(isnan(y)),:) = [];
balldesign(find(isnan(y)),:) = [];
y(find(isnan(y))) = [];
a = (X.'*X)\(X.'*y);
b=a(1);
m=a(2);


[b,bint,r,rint,stats] = regress(y,[num_trls])
[b,bint,r,rint,stats] = regress(y,[ones(numel(like_ratio),1) like_ratio], 0.05)

% Diagnose outliers by finding the residual intervals rint that do not contain 0.

contain0 = (rint(:,1)<0 & rint(:,2)>0);
idx = find(contain0==false)

hold on
scatter(y,r)
scatter(y(idx),r(idx),'b','filled')
xlabel("accuracy")
ylabel("Residuals")
hold off


boxplot(y, balldesign)
boxplot(y, subject)
ylabel('perf 30-150ms')
boxplot(y, pattern)
boxplot(y, farbe)
boxplot(y, like_dich)
ylabel('perf 30-150ms')
xlabel('1 = like, 2 = dislike, 3 = dontcare')
y_like_dislike=y;

ind_3 = find(like_dich == 3);

like_dich(ind_3)=[];
y_like_dislike(ind_3)=[];
y_like_dislike(find(isnan(y_like_dislike))) = [];
like_dich(find(isnan(y_like_dislike))) = [];

[p,tbl,stats] = anovan(y_like_dislike, like_dich );



ind_like = find(like_dich==1);
ind_dislike = find(like_dich==2);
perf_like = y_like_dislike(ind_like);
perf_dislike = y_like_dislike(ind_dislike);

figure;
mean_perf_like = mean(perf_like)
mean_perf_dislike = nanmean(perf_dislike)

figure
scatter(y1,y2)
% non-linear regression:

[beta,R,J,CovB,MSE,ErrorModelInfo]  = nlinfit(X,Y,modelfun,beta0)





end


%% http://compneurosci.com/wiki/images/e/e6/Repeated_ANOVA_MATLAB_v2.pdf   
Biases = rand(18,3*5*2); % subjects, HRs*obstacle pos.*VFs
varNames = cell(3*5*2,1);
for i = 1 : 3*5*2
 v = strcat('V',num2str(i));
 varNames{i,1} = v;
end
% Create a table storing the respones:
tbiases = array2table(Biases, 'VariableNames',varNames);    
Create a table reflecting the within subject factors
HRs = cell(3*5*2,1); % head roll conditions
VFs = cell(3*5*2,1); % Visual feedback conditions
OPs = cell(3*5*2,1); % Obstacle Positions    
    
% Assiging the values to the parameters based on the data sorting
c1 = cell(1,1); c1{1} = 'Y'; c1 = repmat(c1,15,1); VFs(1: 15,1) = c1;
c1 = cell(1,1); c1{1} = 'N'; c1 = repmat(c1,15,1); VFs(16: end,1) = c1;

c1 = cell(1,1); c1{1} = 'HR0'; c1 = repmat(c1,10,1); HRs(1:3:end,1) = c1;
c1 = cell(1,1); c1{1} = 'HRL'; c1 = repmat(c1,10,1); HRs(2:3:end,1) = c1;
c1 = cell(1,1); c1{1} = 'HRR'; c1 = repmat(c1,10,1); HRs(3:3:end,1) = c1;

for i = 1 : 5
 o = strcat('O',num2str(i));
 c1 = cell(1,1); c1{1} = o; c1 = repmat(c1,3,1); OPs((i-1)*3+1:i*3,1) = c1; 
end

OPs(16:end,1) = OPs(1:15,1);

% Create the within table
factorNames = {'HRs','VisualFeedback', 'ObstaclePos'};
within = table(HRs, VFs, OPs, 'VariableNames', factorNames);
% fit the repeated measures model
rm = fitrm(tbiases,'V1-V30~1','WithinDesign',within);
[ranovatblb] = ranova(rm, 'WithinModel','HRs*VisualFeedback*ObstaclePos');
