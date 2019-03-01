% Sourcedummy in this case is a source structure containing the grid .pos
% and .inside field

grad = [];
grad.gradpos = chanpos;
grad.chanpos = chanpos;
tmp = num2cell(chanpos);
% grad.label = label';
for o=1:length(atlas_cluster)
    label_short{o,1} = num2str(atlas_cluster{o,1});
end
grad.label = label_short;

cfg = [];
cfg.grad = grad;
cfg.projection = 'polar';
layout = ft_prepare_layout(cfg);

cfg=[];
cfg.layout = layout;
ft_layoutplot(cfg);