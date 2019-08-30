function avg = fte_subaverage(cfg, data)
    
    if ~isfield(cfg, 'mode') || strcmp(cfg.mode, 'consecutive')
        avg = consecutive_average(cfg, data);
    elseif isfield(cfg, 'mode') && strcmp(cfg.mode, 'bootstrap')
        avg = bootstrap_average(cfg, data);    
    end


end

function avg = consecutive_average(cfg, data)
    
    if isfield(cfg, 'trials')
        inds = cfg.trials;
    else
        inds = 1:numel(data.trial);
    end

    num = length(inds) / cfg.averages;
    trials = cell(1, num);
    for n=1:num
        tmpcfg = [];
        sel = (((n-1)*cfg.averages)+1):(n*cfg.averages);
        tmpcfg.trials = inds(sel);
        tmp = ft_timelockanalysis(tmpcfg, data);
        trials{n} = tmp.avg;
    end

    avg = tmp;
    avg.trial = trials;
    avg = rmfield(avg, {'avg', 'dof', 'var'});
    time = avg.time;
    avg.time = cell(1, size(avg.trial, 2));
    for n=1:size(avg.trial, 2)
        avg.time{n} = time;
    end
end

function avg = bootstrap_average(cfg, data)
    
    if isfield(cfg, 'trials')
        inds = cfg.trials;
    else
        inds = 1:numel(data.trial);
    end
    

    
    trials = cell(1, cfg.repetitions);
    total = length(inds);
    for n=1:cfg.repetitions
        tmpcfg = [];
        sel = randi(total, cfg.averages, 1);
        tmpcfg.trials = inds(sel);
        tmp = ft_timelockanalysis(tmpcfg, data);
        trials{n} = tmp.avg;
    end
    
    avg = tmp;
    avg.trial = trials;
    avg = rmfield(avg, {'avg', 'dof', 'var'});
    time = avg.time;
    avg.time = cell(1, size(avg.trial, 2));
    for n=1:size(avg.trial, 2)
        avg.time{n} = time;
    end
end