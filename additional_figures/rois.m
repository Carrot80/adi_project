
    h(k)=figure;
    scatter3(chanpos(:,1),chanpos(:,2),chanpos(:,3), 'k');
    hold on
   
    fs=fieldnames(indx);
   for k = 1:length(fieldnames(indx))
        h(k)=figure;
        scatter3(chanpos(:,1),chanpos(:,2),chanpos(:,3), 'k');
        hold on
        scatter3(chanpos(indx.(fs{k}),1),chanpos(indx.(fs{k}),2),chanpos(indx.(fs{k}),3), 'g', 'filled');
        h_legend = legend(fs{k}, 'Location','Northeast')  ;  
        hold off
   end
   
   
   savefig(h, [grouppath '_sign_clusters.fig'])
close(h)

   savefig(h, [cd filesep 'rois.fig'])
   
    
    num_rois=unique(cell2mat(atlas_cluster(:,1)));
    p=1;
    for k=1:length(num_rois)
        ind_roi=find(cell2mat(atlas_cluster(:,1))==k);
        scatter3(chanpos(ind_roi,1),chanpos(ind_roi,2),chanpos(ind_roi,3), 'g', 'filled');
        
        
    end
    
    
    for k=1:length(num_rois)
        ind_roi=find(cell2mat(atlas_cluster(:,1))==k);
        scatter3(chanpos(ind_roi,1),chanpos(ind_roi,2),chanpos(ind_roi,3), 'g', 'filled');
        
        
    end
    
    
    
    ind_roi = unique(atl); 
        ind=zeros(1,1);
    for p=1:length(roi)
        [indx.(roi{p})] = find(strcmp(atlas_cluster(:,2), roi{p}));
    end
    
    
    hold on
    scatter3(chanpos(ind,1),chanpos(ind,2),chanpos(ind,3), [], raweffectLike_vs_Dislike.avg(ind, ind_sum_pos_neg(k)), 'filled')
    
    
    
    
    
    
    
    ind = find(pos_neg(:,ind_sum_pos_neg(k)));

    if ~isempty(roi)
        roi_num = cell2mat(atlas_cluster(vox_index(ind,1)));           
    else
        roi_num = cell2mat(atlas_cluster(ind,1));
    end

    if 1 == length(unique(roi_num))
        temp = find(cell2mat(atlas_cluster(:,1))== unique(roi_num));
        roi_name = atlas_cluster(temp(1),2);
        scatter3(chanpos(ind,1),chanpos(ind,2),chanpos(ind,3), [], raweffectLike_vs_Dislike.avg(ind, ind_sum_pos_neg(k)), 'filled')
        h_legend = legend({'brain', char(roi_name)}, 'Location','Northeast')  ;  
        clearvars roi_name  
    else 
        for i=1:length(unique(roi_num))
            unique_roi = unique(roi_num);
            temp = find(cell2mat(atlas_cluster(:,1))== unique_roi(i));
            roi_name(i) = atlas_cluster(temp(1),2);
            clearvars temp
        end

        scatter3(chanpos(ind,1),chanpos(ind,2),chanpos(ind,3), [], raweffectLike_vs_Dislike.avg(ind, ind_sum_pos_neg(k)), 'filled')
        h_legend = legend({'brain', char(roi_name)}, 'Location','Northeast') ;   
        clearvars roi_name  
    end

    title([balldesign ' ' num2str(sign_time(k)) ' s' ])
%     set(h_legend,'Position',[0 0 0.4429 0.3349])
    drawnow
    ax = gca;
    ax.Units = 'pixels';
    ax.SortMethod ='ChildOrder'; 
    position = ax.Position;
    ti = ax.TightInset;
    rect = [-ti(1), -ti(2), position(3)+ti(1)+ti(3), position(4)+ti(2)+ti(4)];
    M(k)=getframe(ax,rect);

    if k==1 
        export_fig ([ grouppath filesep balldesign '_sign_clusters.pdf']);
    else
        export_fig ([grouppath filesep balldesign '_sign_clusters.pdf'], '-append')
    end