function [hypimg]=ImageReconstruction(Sigma1,Y1,meanst,P,tf,params)

[Wrm,e1]=svd(Sigma1);
%e1=real(e1);
Wrm=real(Wrm);
W2=Wrm(:,1:5);%rank(Sigma1));
Y_est = cell(1,params.patterns);
for i=1:params.patterns
    Y_est{i}=pinv(P{i}'*W2)*Y1{i};%
end
Y_est=W2*cell2mat(Y_est);%
Ys1 = zeros(floor(params.L/params.disminuirbands), params.M/params.disminuir*params.N/params.disminuir);
n = 1;
for i=1:params.patterns
    Ys1(:,tf==i)=(Y_est(:,n:((n-1)+sum(tf==i))));
    n = n + sum(tf==i);
end

Ys1=Ys1+meanst;
hypimg = reshape(Ys1',[params.M/params.disminuir,params.N/params.disminuir,floor(params.L/params.disminuirbands)]);

