function formatSensingMatrix()
%% Format Sensing matrix
% Note that the measurements are given by y=Hf=diag(P_i)f, hence, the matrix H
% must be characterized in the lab. But, the camera record the whole H
% hence, this algorithm compare all the generated filters in order to group
% the similar ones into a single P_i as the mean of the filters. This
% should be run only once if the code aperture does no vary

% we capture 20 code apertures, but we will use only 8 out of the 20.
nshots = [9,13,10,14,11,15,12,16];

% the size of the sensor is 776x1032 but we take a smaller window
M = 712;
N = 1000;
L = 101;

% starting pixel of the selected window
start1 = 16;
start2 = 25;

% We designed the coded paerture tu generate 32 partitions, but given the 
% difference in the pixel pitch between the coded aperture and the sensor
% it could be different, hence we are looking for 64.
patterns = 32;

% A decimation factor (Due to the micromirror device artifacts, it is
% better to deciamte the measurements)
disminuir = 2;

% contains the illumination distribution given by the monocromator, it will
% used to remove the difference of the illumination in the coded apertures.
load('codedAperture/meanX')

shots=length(nshots); 


% constant for decimation in the spectral axis, we measure 101 spectral
% bands, but the prism only produce 37
disminuirbands=2.7;

% Preallocate
H = single(zeros(M/disminuir,N/disminuir,floor(L/disminuirbands),shots));
for k=1:shots
    % load a coded aperture
    load(['codedAperture/',num2str(nshots(k)),'_aperturefinalfinal_raro'])
    
    % divide by the illumination (to remove the illumination difference)
    X = X ./ meanX;

    % select the desired window of the measurements.
    sH = X(start1:start1+M-1,start2-14:start2-14+N-1,:);%
    
    % preallocate to store the decimated version
    G=zeros(round(M/disminuir),(N/disminuir),L);
    for si=1:L
        G(:,:,si) = imresize(sH(:,:,si),[round(M/disminuir),(N/disminuir)],'bilinear');
    end

    % We measure the spectral response of the coded aperture every 2 nm, however, due
    % to the dispersion of the prism, we reduce the number of the bands
    % following this dispersion (which is not linear)
    if disminuirbands>1
        [bandsdisc,G1] = selectBands(G,disminuirbands,M,N,disminuir);
    else
        G1=G;
    end
    % Increase the contrast on the coded aperture
    G1 = G1.^1.1;
    % normalize to 1
    G1=G1./max(max(G1,[],2),[],1);
    H(:,:,:,k)=single(G1);
end
clear sH X G1 G

[tf,distances]=locateFilterDistribution(H,patterns,M,N,L,disminuir,shots,disminuirbands);
% place each shot in a cell
Hs = cell(1,shots);
for k=1:shots
    % place each shot of the coded aperture in the cell in matrix form
    Hs{k} = ((reshape(H(:,:,:,k),[M/disminuir*N/disminuir,floor(L/disminuirbands)]))')';
end

% create the sensing matrices P
P = cell(1,patterns);
for i=1:patterns
    P{i}=zeros(shots,floor(L/disminuirbands));
    for j=1:shots
        % using the spatial distribution (tf) to select the similar filters
        % and select the median as the sensing filter for all this pixels
        % in the image
        P{i}(j,:) = median(Hs{j}(tf==i,:),1);
    end
end
%save('ready/H','H')
save('ready/tf','tf')
save('ready/P','P')
save('ready/nshots','nshots')
params.M=M;
params.N=N;
params.L=L;
params.start1=start1;
params.start2=start2;
params.disminuir=disminuir;
params.disminuirbands=disminuirbands;
params.patterns=patterns;
save('ready/params','params')

%save('ready/df','df')