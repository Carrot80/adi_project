
function [var, MEG_interp] = checkfields (var, MEG_interp)

    if ~isfield(MEG_interp, 'dimord')
        MEG_interp = setfield(MEG_interp, 'dimord', 'chan_time');
    end
    if ~isfield(var, 'dimord')
        var = setfield(var, 'dimord', 'chan_time');
    end
    if ~isfield(MEG_interp, 'sampleinfo')
        MEG_interp.sampleinfo = MEG_interp.sampleinfo_orig;
    end
    if ~isfield(var, 'sampleinfo')
        var.sampleinfo = var.sampleinfo_orig;
    end

    fieldnames1 = fieldnames(var);
    fieldnames2 = fieldnames(MEG_interp);
    [field, row] = setdiff(fieldnames1, fieldnames2);
    var = rmfield(var, field);
    clear field row fieldnames1 fieldnames2
    fieldnames1 = fieldnames(var);
    fieldnames2 = fieldnames(MEG_interp);
    [field, row] = setdiff(fieldnames2, fieldnames1);
    MEG_interp = rmfield(MEG_interp, field);
    MEG_interp = orderfields(MEG_interp, var);
    
end