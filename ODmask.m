function [gr_img1 a b]=ODmask(img)
res_img = imresize(img,[256 256]);
gr_img1= res_img;
h = fspecial('disk',5);
filt_img = imfilter((res_img-40),h);
bw = im2bw(filt_img,0.45);
se = strel('disk',5);
BW_erod = imopen(bw,se);
se = strel('disk',9);
BW_dil = imdilate(BW_erod,se);
P=bwlabel(BW_erod,8);
BB=regionprops(P,'Boundingbox');
BB1=struct2cell(BB);
BB2=cell2mat(BB1);
[s1 s2]=size(BB2);
mx=0;
for k=3:4:s2-1
    p=BB2(1,k)*BB2(1,k+1);
    if p>mx & (BB2(1,k)/BB2(1,k+1))<1.8
        mx=p;
        o=k;
    end
end
i1=round(BB2(1,o-2))-1;
i2=round(BB2(1,o-1))-2;
i3=BB2(1,o)+10;
i4=BB2(1,o+1)+10;
if i3<30
    i1=i1-7;
    i2=i2-7;
    i3=i3+10;
    i4=i4+10;
end
a=[i1 i2 i3 i4];
i1=round(BB2(1,o-2))-20;
i2=round(BB2(1,o-1))-20;
i3=BB2(1,o)+40;
i4=BB2(1,o+1)+40;
if i3<70
    i1=i1-20;
    i2=i2-20;
    i3=i3+40;
    i4=i4+40;
end

if i4>60
    i4=i4-20;
end
b=[i1 i2 i3 i4];