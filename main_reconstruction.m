if isfile('ready/P.mat') && isfile('ready/tf.mat')
     load('ready/P.mat')
     load('ready/tf.mat')
     load('ready/params')
     load('ready/nshots.mat')
else
     error('P and tf matrices does not exist, please run formatSensingMatrix.m first')
     return;
end

% check that the sensing matrix dimensions match
[s1,s2]=size(P{1});
if (s2>s1)
    for i=1:params.patterns
        P{i}=P{i}';
    end
end
P1=P;

nshots1=nshots;
% got better result if I reconstruct 2 shots eachtime and average
% reconstruction
shots = 2;
meanst1=0;
shotstotal=length(nshots1);

% store average covariance
Sigmatotal=0;
hypimgtotal=0;
for i=1:(shotstotal/shots)
    nshots=nshots1((1:shots)+((i-1)*shots))
    for r=1:params.patterns
        P{r}=P1{r}(:,(1:shots)+((i-1)*shots));
    end
    load('measurements/paisa_aperturefinalfinal.mat')
    for jjk=1:16
        X(:,:,jjk)=(X(:,:,jjk)./2.55);
    end
    X=imrotate(X,-0.25,'bilinear');
    % select the same region that in the formatSensingMatrix script
    sH1 = X(params.start1:params.start1+params.M-1,params.start2:params.start2+params.N-1,:);
    % reshape the measurements to a matrix form
    F=double(reshape(sH1(:,:,nshots),[params.M*params.N*shots,1]));
    Y = zeros(round(params.M/params.disminuir),(params.N/params.disminuir),shots);
    for si=1:shots
        % decimate the measurements as the coded aperture
        Y(:,:,si) = imresize(sH1(:,:,nshots(si)),[round(params.M/params.disminuir),(params.N/params.disminuir)],'bilinear');
    end
    % reshape to a matrix form
    Y1 = reshape(Y,[round(params.M/params.disminuir)*round(params.N/params.disminuir),shots]);
    % Spatial distribution of the coded aperture
    tf = reshape(tf,[round(params.M/params.disminuir)*round(params.N/params.disminuir),1]);
    Ys = cell(1,params.patterns);
    for i=1:params.patterns
        % create the subsets (partitions)
        Ys{i}=Y1(tf==i,:)';
    end

    

    meanst=estimate_mean(P,Ys);
    Y1=cell(size(Ys));

    for r=1:params.patterns
        Y1{r}=Ys{r}-P{r}'*(kron(meanst,ones(1,size(Ys{r},2))));%meanst
    end

    for j=1:15
        mst=estimate_mean(P,Y1);
        %mst = max(mst,0);
        meanst=meanst+mst;
        for r=1:params.patterns
            Y1{r}=Y1{r}-P{r}'*(kron(mst,ones(1,size(Ys{r},2))));
        end
    end
    meanst1 = meanst1+meanst;
    Sigma1 = estimate_cov6(P,Y1,240,0.0002,200,1.5,0,0,0);
    Sigmatotal = Sigmatotal + Sigma1;
    hypimg=ImageReconstruction(Sigma1,Y1,meanst,P,tf,params);
    hypimgtotal = hypimgtotal + hypimg;
end
meanst1=meanst1./(shotstotal/shots);
Sigmatotal=Sigmatotal./(shotstotal/shots);
hypimgtotal = hypimgtotal/(shotstotal/shots);
imrgb = RGBcomposite(hypimgtotal,params);

imagesc(imrgb)
