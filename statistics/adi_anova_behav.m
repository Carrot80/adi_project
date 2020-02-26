for kk=1:length(balldesign)
    if contains(balldesign{kk}, 'gb')
        color{kk,1} = 'gb';
    elseif contains(balldesign{kk}, 'gg')
        color{kk,1} = 'gg';
    elseif contains(balldesign{kk}, 'rw')
        color{kk,1} = 'rw';
    end
    if contains(balldesign{kk}, 'f')
        pattern{kk,1} = 'foot';
    elseif contains(balldesign{kk}, 's')
        pattern{kk,1} = 'space';
    elseif contains(balldesign{kk}, 'v')
         pattern{kk,1} = 'volley';
    end
end

stat = anova()
p = anovan(like,{color,pattern})
[p,tbl,stats,terms] = anovan(like,{pattern, color},'model','interaction','varnames',{'pattern','color'})
[p,tbl,stats,terms] = anovan(like,{color, pattern},'model','interaction','varnames',{'color', 'pattern'})

boxplot(like, pattern)
boxplot(like, color)

like_foot = like(find(strcmp(pattern, 'foot')));
like_volley = like(find(strcmp(pattern, 'volley')));
like_space = like(find(strcmp(pattern, 'space')));
[h,p,ksstat,cv]  = kstest(like_space)
[c,m,h,s] = multcompare(stats,'Alpha',0.05,'CType','bonferroni','Display','on')
[c,m,h,s] = multcompare(stats,'Alpha',0.05,'CType','bonferroni', 'Estimate','row')
