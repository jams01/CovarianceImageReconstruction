function [Sigma1,Yn, var_noise, minerror] = estimate_cov6(P,Y,it,tau,noise,softp,D1,reg,Sigmareal)
%% Estimate covariance matrix
% Input:
% P: is a cell array with sensing matrices such as Y{i}=P{i}'*F{i};
% Y: is a cell array with measurements.
% it: is the number of iterations
% nev: is the rank of the covariance matrix.
% noise: SNR of the measurements.
% Sigmareal: is the original covariance matrix (no used by algorithm, only
% for convergence analysis).
% soft: use a filter step in the gradient.
% reg: regularization type : 0 -> trace (low rank), 1 -> toeplitz 
%
% Output:
% Sigma1: estimated covariance matrix.
% Yn : signal with noise
%%
partition=length(Y);
D=cell(1,partition);
%k=zeros(1,partition);
%Sigmar=0;
Sigmas=cell(1,partition);
Yn=cell(1,partition);
var_noise = 0;
for i=1:partition %realiza el proceso de muestreo y usa los valores singulares para mejor inversion
%     y0=Y{i};
%     Y{i} = awgn(Y{i},noise,'measured'); %agrego ruido
%     Yn{i} = Y{i};
%     temp = Yn{i} - y0;
%     var_noise = var_noise + trace(temp*temp'./size(temp,2))/size(y0,1);
    Sigmas{i}=Y{i}*Y{i}'./size(Y{i},2);
end
var_noise = var_noise / partition;
Y1=Y;
Di=cell(1,length(D));
Dj=cell(1,length(D));
Pj=cell(1,length(D));

for i=1:partition
    Di{i}=pinv(P{i}*P{i}')'*P{i};
    Dj{i}=pinv(P{i}'*P{i})*P{i}';
    Pj{i}=pinv(P{i});
end
partition=length(D);

Sigma1=0;
%grad1=cell(partition,1);
for i=1:1:partition
    %grad1{i} = ((((pinv(P{i}')*Y1{i}))*(Y1{i}'*pinv(P{i}')'))/size(Y1{i},2));
    Sigma1=Sigma1+((((pinv(P{i}')*Y1{i}))*(Y1{i}'*pinv(P{i}')'))/size(Y1{i},2)); %inicialiso mi sigma
    %grad1=grad1+(P{(i)}*P{(i)}');
end

Sigma1=Sigma1./(partition*0.5);


wnd = 4;
nold1=zeros(1,wnd);


Sigmaold=Sigma1;

I=eye(size(Sigma1));
l=trace(Sigma1)*tau;
nold = calc_obj_fun1(P,Sigma1,Sigmas,1);
nold1(:)=nold;
s=1;
n=4;
L=s;
minerror=100000;
%tau2=0.0;
momentum1 = 0.2;
gradold=0;
for i=1:it
    Sigmao=Sigma1;
    Sigma1=Sigmao;
    tt=randperm(partition);
    grad=zeros(size(Sigma1));

    for j=1:floor(partition/1)
        R=Di{tt(j)}*((P{tt(j)}'*(Sigmao)*P{tt(j)})-Sigmas{tt(j)})*Dj{tt(j)};%+ tau2.* (P{tt(j)}*P{tt(j)}') ;
        Sigmao = Sigmao - ((1/L).*R);%(1/L)
        grad=grad+R;
    end
    L=s;
    %grad=grad./(64*64);
    %grad=grad-D1;
    if(softp>0)
        grad=imgaussfilt(grad-D1,softp);
    else
        grad=grad-D1;
    end
    Sigmaold = Sigma1;
    it1=8; % OJO ACA
    if(reg==0)
        for k=0:it1
            z=Sigma1 - (1/L)*(((1-momentum1).*grad + momentum1.*gradold) + l*I);%1/L
            [W,e] = eig((z)./1);
            e=diag(max(diag((e)),0));
            z=real(W*e*W');
            obj_f_l = calc_obj_fun1(P,z,Sigmas,l);
            obj_x_k = calc_obj_fun1(P,Sigma1,Sigmas,l);
            ob2 = trace((grad+ l*I)'*(z-Sigma1));
            %ob2 = sum(sum((grad+ l*I).*(z-Sigma1)));
            %ob3 = L/2*sum(sum((z-Sigma1).*(z-Sigma1)));
            ob3 = L/2 *norm((z-Sigma1),'fro')^2;
            if(obj_f_l<=1.01*(obj_x_k + ob2 + ob3))
                %l=trace(Sigma1)*tau;
                Sigma1=z;
                break;
                
            elseif(k==it1)
                Sigma1=z;
            else
                L=L*n;
            end
        end

    end
    s=L;
    %l=trace(Sigma1)*tau;
    
    gradold=grad;
    n2 = calc_obj_fun1(P,Sigma1,Sigmas,1);
    nold1=[nold1(2:end),n2];
    if (norm(Sigma1-Sigmareal,'fro')<minerror)
       minerror =  norm(Sigma1-Sigmareal,'fro');
    end
    if (mod(i,wnd)==0 || i==1)
        if Sigmareal==0
            fprintf('it=%d est=%f sigmas=%d  rank=%d, L=%d\n',i,n2,norm(Sigma1-Sigmaold,'fro')./norm(Sigma1,'fro'),rank(Sigma1),L)
        else
            fprintf('it=%d est=%f sigmas=%d  rank=%d, L=%d, real=%.6f\n',i,n2,norm(Sigma1-Sigmaold,'fro')./norm(Sigma1,'fro'),rank(Sigma1),L, norm(Sigmareal-Sigma1,'fro'))
        end
    end
    if(norm(Sigma1-Sigmaold,'fro')/norm(Sigma1,'fro')<2e-9)
        fprintf('it=%d est=%f sigmas=%d  rank=%d, L=%d\n',i,n2,norm(Sigma1-Sigmaold,'fro')./norm(Sigma1,'fro'),rank(Sigma1),L)
       break; 
    end
end
clear Sigmas
clear Di
clear Dj
clear Pj
end