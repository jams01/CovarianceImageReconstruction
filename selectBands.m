function [bandsdisc, G1]=selectBands(G,disminuirbands,M,N,disminuir)
% Preallocate
G1 = zeros(M/disminuir,N/disminuir,floor(101/disminuirbands));
index=1;
% This are the captured bands in the lab
bandsall=450:2:650;
% The selected bands
bandsdisc = zeros(1,floor(101/disminuirbands));
% Measured dispersion in the lab
load('disper')
for i=1:floor(101/disminuirbands)
    G1(:,:,index) = sum(G(:,:,disper==i),3)./sum(disper==i);
    [v,p]=find(disper==i);
    bandsdisc(i)=mean(bandsall(p));
    index=index+1;
end