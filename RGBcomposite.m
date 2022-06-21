function [imt]=RGBcomposite(imhyp,params)

%simulate a integration response of a camera
integration = normpdf([-2:.20:2.1],0,1);
integration=integration./sum(integration);
integration(:) = integration / sum(integration);
% we recover the spectral signature of white and store it to normalize the
% data
load('b41')
b=reshape(b,[1,1,floor(params.L/params.disminuirbands)]);
%normalize the data w.r.t the white (white balance)
imhyp1=imhyp;
for i=1:params.M/params.disminuir
    for j=1:params.N/params.disminuir
        imhyp1(i,j,: )=imhyp(i,j,:)./(b);
    end
end
imt =zeros(params.M/params.disminuir,params.N/params.disminuir,3);
% integrate
int1 = integration(3:20)./sum(integration(3:20));
int1(1:3)=0;
imt(:,:,3) =  sum(imhyp1(:,:,1:18).*reshape(int1,[1,1,18]),3);
int1 = integration(3:17)./sum(integration(3:17));
imt(:,:,2) =  sum(imhyp1(:,:,17:31).*reshape(int1,[1,1,15]),3);
int1 = integration(5:13)./sum(integration(5:13));
imt(:,:,1) =  sum(imhyp1(:,:,29:end).*reshape(int1,[1,1,9]),3);

imt = imt/1.5;