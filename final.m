clc
close all
clear all
[f p]=uigetfile('*.jpg','*.tiff','Select An Image');
I = imread([p f]);
I=imresize(I,[256 256]);
m = zeros(size(I,1),size(I,2));
G=I(:,:,2);
R=I(:,:,1);
ig=G;
[r c]=size(G);
for i=1:r
    for j=1:c
        if G(i,j)>5
            G(i,j)=255;
        end
    end
end
figure(1);subplot(3,3,1);
imshow(G);axis off;
title('Image background');
figure(1);subplot(3,3,2);
imshow(ig);axis off;
title('Green Channel Image');
 h=ig;
 k=ig-100;
 k=imresize(k,[256 256]);
 Threshold = 10;
 b = boundary(k, Threshold);
figure(1);subplot(3,3,3);
imshow(b);axis off;
title('OD Region');
img=imresize(ig,[256 256]);
imS = medfilt2(img, [10 10]);
imS=imS-10;
figure(1);subplot(3,3,4);imshow(imS);axis off;title('Optical cup Segmentation');
[OD a b1]=ODmask(img);
i1=a(1);
i2=a(2);
i3=a(3);
i4=a(4);
figure(1);subplot(3,3,5);imshow(OD);axis off;
m(i2:(i2+i4),i1:(i1+i3))=1;
imshow(m);axis off;
 seg = odsegment(I,m,500,0.1,'chan');
 i1=b1(1);
 i2=b1(2);
 i3=b1(3);
 i4=b1(4);
 ig=imresize(R,[256 256]);
i5=i1+i3;
i6=i2+i4;
if i5>=size(R,1)
    i5=size(R,1);
end
[Im0] =ig(i2:i6,i1:i5);
Im0=imresize(Im0,[100 100]);
[Ny,Nx] = size(Im0);
Im0 = double(Im0);
Im0 = ( Im0-min(Im0(:)) )/ ( max(Im0(:))-min(Im0(:)) );
c1 = 0.26;
c2 = 0.91;
wrin = (Im0 - c2).^2;
wrout = (Im0 - c1).^2;
wr = wrout - wrin;
Gauss = fspecial('gaussian',[Nx Ny],2);
Im0s = real(ifftshift(ifft2( fft2(Im0).* fft2(Gauss) )));
[Gx,Gy] = gradient(Im0s);
NormGrad = sqrt(Gx.^2 + Gy.^2); 
NormGrad = sqrt(Gx.^2 + Gy.^2); 
wb = 1./ (1 + 30* abs(NormGrad).^2); 
phi0 = ones(size(Im0));
w=9; phi0(1:end,1:end)= -1;
phi0 = FMReDistSDF_mex(single(phi0),single(sqrt(Nx^4+Ny^4)));
phi0 = double(phi0);
figure(1);subplot(3,3,6);
axis off;
imagesc(Im0); colormap('gray'); title('Image');axis off;
phi  = Segment2D_public(wb,wr,-4,phi0,60);
figure(1);subplot(3,3,7);
axis off;
imagesc(Im0); colormap('gray'); hold on; contour(phi,[0 0],'m'); title('final segmentation');
for i=1:100
for j=1:100
if phi(i,j)<0
iot(i,j)=1;
else
iot(i,j)=0;
end
end
end
k1=1;
 for i=1:100
for j=1:100
if iot(i,j)==1
k1=k1+1;
end
end
end
 ig=imresize((imS+20),[256 256]);
i5=i1+i3;
i6=i2+i4;
if i5>=size(R,1)
    i5=size(R,1);
end
[Im0] =ig(i2:i6,i1:i5);
Im0=imresize(Im0,[100 100]);
[Ny,Nx] = size(Im0);
Im0 = double(Im0);
Im0 = ( Im0-min(Im0(:)) )/ ( max(Im0(:))-min(Im0(:)) );
c1 = 0.26;
c2 = 0.91;
wrin = (Im0 - c2).^2;
wrout = (Im0 - c1).^2;
wr = wrout - wrin;
Gauss = fspecial('gaussian',[Nx Ny],2);
Im0s = real(ifftshift(ifft2( fft2(Im0).* fft2(Gauss) )));
[Gx,Gy] = gradient(Im0s);
NormGrad = sqrt(Gx.^2 + Gy.^2); 
% NormGrad = sqrt(Gx.^2 + Gy.^2); 
wb = 1./ (1 + 30* abs(NormGrad).^2); 
phi0 = ones(size(Im0));
w=9; phi0(1:end,1:end)= -1;
phi0 = FMReDistSDF_mex(single(phi0),single(sqrt(Nx^4+Ny^4)));
phi0 = double(phi0);
phi  = Segment2D_public(wb,wr,-4,phi0,60);
figure(1);subplot(3,3,7); hold on; contour(phi,[0 0],'v'); title('Final segmentation');
axis off;
for i=1:100
for j=1:100
if phi(i,j)<0
iot(i,j)=1;
else
iot(i,j)=0;
end
end
end
k2=1;
 for i=1:100
for j=1:100
if iot(i,j)==1
k2=k2+1;
end
end
end
 CuptoDiscRatio=k2/k1/2
 if CuptoDiscRatio<0.3
     msgbox('Normal image');
 else
     msgbox('Glaucoma image');
 end
R=I(:,:,1);
G=I(:,:,2);
B=I(:,:,3);
a=imhist(R);
b=imhist(G);
c=imhist(B);
Testfeature(1,1:256)=a; 
Testfeature(1,257:512)=b; 
Testfeature(1,513:768)=c; 
save Testfeature Testfeature
for i=1:20
path = 'Dataset\';
ext = '.jpg';
file = strcat(path,num2str(i),ext);
img=imread(file);
img=imresize(img,[256 256]);
R=img(:,:,1);
G=img(:,:,2);
B=img(:,:,3);
at(i,:)=imhist(R);
bt(i,:)=imhist(G);
ct(i,:)=imhist(B);
Trainfeature(i,1:256)=at(i,:); 
Trainfeature(i,257:512)=bt(i,:); 
Trainfeature(i,513:768)=ct(i,:); 
save Trainfeature Trainfeature
end
load Truelabel
class=svm(round(Testfeature));
if class == 1
msgbox('SVM Classifies as Normal Image');
elseif class == 2
msgbox('SVM Classifies as Glaucoma Image');
end
tclass = svm(Trainfeature);
