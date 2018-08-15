function pred=gnbpredict(y, x, model)
    [M,N]=size(x);

    pred = nan(M, 1);
    for i=1:M
        lhy1=model.prior(1);
        lhy0=model.prior(2);
        for j=1:N
            lhy1=lhy1/sqrt(model.cvar(1,j))*exp(-(x(i,j)-model.cmean(1,j))*(x(i,j)-model.cmean(1,j))/2/model.cvar(1,j));
            lhy0=lhy0/sqrt(model.cvar(2,j))*exp(-(x(i,j)-model.cmean(2,j))*(x(i,j)-model.cmean(2,j))/2/model.cvar(2,j));
        end
        if lhy1>lhy0
            pred(i)=model.classes(2);
        else
            pred(i)=model.classes(1);
        end
    end
end