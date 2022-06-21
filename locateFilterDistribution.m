function [tf,distances]=locateFilterDistribution(H, patterns, M, N, L, disminuir,shots,disminuirbands)
%% Given the coded aperture, this function find the position where the coded pattern is equal
% i.e. the coded patters are repetead in the image but we do not know what codes are they and where they
% are
% input: 
%      H: coded aperture 
%      patters: number of patterns that whant to look for
% output:
%      tf: spatian distribution of the codes in H
%      distances: distance measure of each pixel to the reference patterns

distances1=0;
row = 100;
for j=1:shots
    % select the reference patterns, theoretically, each row has the same
    % filters, but we select the 100th row
    t1 = [squeeze(H(100,24:(24+patterns),:,rrr))];%H(70,19:19+patterns-1,:,3);%3
    % convert the coded aperture into a matrix, where the columns are the
    % spectral patterns
    H1=reshape(H(:,:,:,j),[M/disminuir*N/disminuir,floor(L/disminuirbands)]);%3 (nada4)
    H1 = normc(H1')';
    distances = zeros(M/disminuir*N/disminuir,patterns);
    for i=1:patterns
        t = t1(i,:);
        t = t ./norm(t);
        % compute the distance between the reference pattern and each pixel
        % in the coded aperture (if the dot product is 1, it means that the
        % patterns are the same)
        distances(:,i) = H1*t';
    end
    distances1=distances1+distances;
end
distances = distances1/shots;
distances = distances + (distances>0.2);
distances =[ones(M*N/disminuir^2,1),distances];
[df,tf] = max(distances');
% spatial distribution of the patterns
tf=tf-1;