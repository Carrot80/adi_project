

%% euklidscher Abstand:

p = [5.22	6.99	1.3];
q = [5.25	6.97	1.33];

euclidsche_norm = sqrt( (q(1)-p(1))^2 + (q(2)-p(2))^2 + (q(3)-p(3))^2);


for ii=1:size(q,1)
    euclidsche_norm(ii,1) = sqrt( (q(ii,1)-p(ii,1)).^2 + (q(ii,2)-p(ii,2)).^2 + (q(ii,3)-p(ii,3)).^2);
end



