
function [] = adi_figureTPRcrossval (Condition1vs2, time, cond1, cond2, NameCond1vs2, outPath, freqbandname, run)
 
figure
plot(time, Condition1vs2.Accuracy);
hold on

indSig = []; 
indSig(1:length(Condition1vs2.Binominal)) = NaN;
sig = find(Condition1vs2.Binominal <= 0.1);
indSig(sig) = Condition1vs2.Binominal(sig);
plot(time, indSig,'r+'); 

title(strcat(NameCond1vs2, '_', freqbandname));
xlabel('time');
ylabel('accuracy/p-value'); 
ylim([0 1]);

savefig([outPath, NameCond1vs2, run, '_', freqbandname,'.fig'])
fig = ([outPath, NameCond1vs2, run, '_', freqbandname]);
print('-dpng', fig); 

total_cond1 = Condition1vs2.Contingency{1,1}(1,1) + Condition1vs2.Contingency{1,1}(1,2);
total_cond2 = Condition1vs2.Contingency{1,1}(2,1) + Condition1vs2.Contingency{1,1}(2,2);

for i = 1:length(Condition1vs2.Contingency)
    tpr_cond1(i)= Condition1vs2.Contingency{1,i}(1,1)/total_cond1; % true positive rate
    tpr_cond2(i)= Condition1vs2.Contingency{1,i}(2,2)/total_cond2; % true positive rate
end

figure
plot(time, Condition1vs2.Accuracy, 'r');
hold on
plot(time, indSig,'r+');  
hold on
plot(time, tpr_cond1,'b'); 
hold on
plot(time, tpr_cond2,'k'); 
leg_cond1 = (['TPR_', cond1]);
leg_cond2 = (['TPR_', cond2]);
legend ({'total TPR', 'p-value of total TPR', leg_cond1, leg_cond2});
title([NameCond1vs2, '_', freqbandname]);
xlabel('time');
ylabel('accuracy/p-value'); 
ylim ([0 1]);
savefig([outPath, NameCond1vs2, run, '_', freqbandname, '_TPR.fig']);
fig = ([outPath, NameCond1vs2, run, '_', freqbandname, '_TPR']);
print('-dpng', fig); 
close all
    
end