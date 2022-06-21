function obj=calc_obj_fun1(D,Sigma1,Sigmas,lambda)
    
    obj = 0;
    for i=1:length(D)
        obj = obj + norm(D{i}'*Sigma1*D{i}-Sigmas{i},'fro')^2;
    end
    obj=obj./length(D);
    obj=obj+lambda.*trace(Sigma1);
end