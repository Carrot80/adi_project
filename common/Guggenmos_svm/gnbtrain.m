function model = gnbtrain(y, x)
    [M,N]=size(x);
    
    model.classes = unique(y);

    indy0=find(y==model.classes(1));
    indy1=find(y==model.classes(2));
    lenindy0=length(indy0);
    lenindy1=length(indy1);

    model.prior(1)=lenindy1/length(y);   %y=1
    model.prior(2)=lenindy0/length(y);
    model.cmean=zeros(2,N);   %first row: y=1
    model.cvar=zeros(2,N);    %first row: y=1

    model.cmean(1,:)=sum(x(indy1,:))/lenindy1;
    model.cmean(2,:)=sum(x(indy0,:))/lenindy0;

    model.cvar(1,:)=sum((x(indy1,:)-repmat(model.cmean(1,:),lenindy1,1)).*(x(indy1,:)-repmat(model.cmean(1,:),lenindy1,1)))/lenindy1;
    model.cvar(2,:)=sum((x(indy0,:)-repmat(model.cmean(2,:),lenindy0,1)).*(x(indy0,:)-repmat(model.cmean(2,:),lenindy0,1)))/lenindy0;


end