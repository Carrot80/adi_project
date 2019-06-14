function [session]=adi_regroup_session(session)




like.trial=session.trial(find(session.labels==1));
like.time=session.time(find(session.labels==1));
like.response_label = session.response_label(find(session.labels==1));
% like.response = session.response(find(session.labels==1));

dislike.trial=session.trial(find(session.labels==2));
dislike.time=session.time(find(session.labels==2));
dislike.response_label = session.response_label(find(session.labels==2));
% dislike.response = session.response(find(session.labels==2));


session=rmfield(session, 'trial');
session=rmfield(session, 'data');
session=rmfield(session, 'labels');
% session=rmfield(session, 'response');
session=rmfield(session, 'response_label');

session.labels=zeros(1,length(like.trial)+length(dislike.trial));
session.labels(1:length(like.trial))=ones(1,length(like.trial));
session.labels(length(like.trial)+1:end)=2*ones(1,length(dislike.trial));

session.trial=[like.trial dislike.trial];
session.response_label = [like.response_label dislike.response_label];
% session.response = [like.response dislike.response];

for i=1:length(session.trial)
   session.data(i,:,:)=session.trial{i};     
end


clearvars like dislike





end