function mean=estimate_mean(P,Y)
    mean=0;
    number_elements=0;
    for i=1:size(Y,2)
       [S,V,D]=svd(P{i}');
        temp=P{i}*pinv(P{i}'*P{i})*Y{i};
        mean=mean+sum(temp,2);
        number_elements=number_elements+size(Y{i},2);
    end
    mean=mean./(number_elements);
end