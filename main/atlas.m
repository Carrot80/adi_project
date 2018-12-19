    for m = 1:length(atlas.tissuelabel)
        ind = find(cell2mat(atlas.sources_roi_numbers(:,1)) == m);
        ind_atlas.(atlas.tissuelabel{m})=ind; 
    end