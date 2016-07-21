
function seg = odsegment(I,mask,num_iter,mu,method)
  if(~exist('mu','var')) 
    mu=0.2; 
  end
  
  if(~exist('method','var')) 
    method = 'chan'; 
  end


   s = 200./min(size(I,1),size(I,2)); 
   if s<1
       I = imresize(I,s);
   end
  

  if ischar(mask)
      switch lower (mask)
          case 'small'
              mask = maskcircle2(I,'small');
          case 'medium'
              mask = maskcircle2(I,'medium');
          case 'large'
              mask = maskcircle2(I,'large');              
          case 'whole'
              mask = maskcircle2(I,'whole'); 
                   
          case 'whole+small'
              m1 = maskcircle2(I,'whole');
              m2 = maskcircle2(I,'small');
              mask = zeros(size(I,1),size(I,2),2);
              mask(:,:,1) = m1(:,:,1);
              mask(:,:,2) = m2(:,:,2);
          otherwise
              error('unrecognized mask shape name (MASK).');
      end
  else
      if s<1
          mask = imresize(mask,s);
      end
      if size(mask,1)>size(I,1) || size(mask,2)>size(I,2)
          error('dimensions of mask unmathch those of the image.')
      end
      switch lower(method)
          case 'multiphase'
              if  (size(mask,3) == 1)  
                  error('multiphase requires two masks but only gets one.')
              end
      end

  end       

  
switch lower(method)
    case 'chan'
        if size(I,3)== 3
            P = rgb2gray(uint8(I));
            P = double(P);
        elseif size(I,3) == 2
            P = 0.5.*(double(I(:,:,1))+double(I(:,:,2)));
        else
            P = double(I);
        end
        layer = 1;
        
    case 'vector'
        s = 200./min(size(I,1),size(I,2)); 
        I = imresize(I,s);
        mask = imresize(mask,s);
        layer = size(I,3);
        if layer == 1
            display('only one image component for vector image')
        end
        P = double(I);
            
    case 'multiphase'
        layer = size(I,3);
        if size(I,1)*size(I,2)>200^2
            s = 200./min(size(I,1),size(I,2)); 
            I = imresize(I,s);
            mask = imresize(mask,s);
        end
            
        P = double(I);  
    otherwise
        error('!invalid method')
end



switch lower(method)
    case {'chan','vector'}
        
        mask = mask(:,:,1);
        phi0 = bwdist(mask)-bwdist(1-mask)+im2double(mask)-.5; 
       
        force = eps; 
       
          figure();
          imshow(I); title('Input Image');
         contour(flipud(phi0), [0 0], 'r','LineWidth',1); title('initial contour');
           title('Segmentation');
        
          for n=1:num_iter
              inidx = find(phi0>=0);
              outidx = find(phi0<0); 
              force_image = 0;  
              for i=1:layer
                  L = im2double(P(:,:,i)); 
                  c1 = sum(sum(L.*Heaviside(phi0)))/(length(inidx)+eps); 
                  c2 = sum(sum(L.*(1-Heaviside(phi0))))/(length(outidx)+eps);
                  force_image=-(L-c1).^2+(L-c2).^2+force_image; 
                 
              end
      
              force = mu*kappa(phi0)./max(max(abs(kappa(phi0))))+1/layer.*force_image;
             
              force = force./(max(max(abs(force))));
              
              dt=0.5;             
              
              old = phi0;
              phi0 = phi0+dt.*force;
              new = phi0;
              indicator = checkstop(old,new,dt);
              
              if(mod(n,20) == 0) 
                 showphi(I,phi0,n);  
              end;
              if indicator 
                  showphi(I,phi0,n);

                  
                  seg = phi0<=0; 

                   figure();imshow(seg); title('Final OD Segmentation');

                  return;
              end
          end;
          showphi(I,phi0,n);
         
          seg = phi0<=0; 
        figure(); imshow(seg); title('Final OD Segmentation');
    case 'multiphase'
       
        mask1 = mask(:,:,1);
        mask2 = mask(:,:,2);
        phi1=bwdist(mask1)-bwdist(1-mask1)+im2double(mask1)-.5;
        phi2=bwdist(mask2)-bwdist(1-mask2)+im2double(mask2)-.5;       
       
        figure();
        subplot(2,2,1); 
        if layer ~= 1
            imshow(I); title('Input Image');
        else
            imagesc(P); axis image; colormap(gray);title('Input Image');
        end
        subplot(2,2,2);
        hold on
        contour(flipud(mask1),[0,0],'r','LineWidth',2.5); 
        contour(flipud(mask1),[0,0],'x','LineWidth',1);
        contour(flipud(mask2),[0,0],'g','LineWidth',2.5);
        contour(flipud(mask2),[0,0],'x','LineWidth',1);
        title('initial contour');
        hold off
      ; title('Segmentation');        
       
        for n=1:num_iter
             
              nb1 = find(phi1<1.2 & phi1>=-1.2);
              inidx1 = find(phi1>=0); 
              outidx1 = find(phi1<0); 

              nb2 = find(phi2<1.2 & phi2>=-1.2);
              inidx2 = find(phi2>=0);
              outidx2 = find(phi2<0);              

              cc11 = intersect(inidx1,inidx2); 
              cc12 = intersect(inidx1,outidx2); 
              cc21 = intersect(outidx1,inidx2); 
              cc22 = intersect(outidx1,outidx2); 
              
              f_image11 = 0;
              f_image12 = 0;
              f_image21 = 0;
              f_image22 = 0; 
              
              for i=1:layer
                  L = im2double(P(:,:,i)); 
          
              if isempty(cc11)
                  c11 = eps;
              else
                  c11 = mean(L(cc11));
              end
              
              if isempty(cc12)
                  c12 = eps;
              else
                  c12 = mean(L(cc12)); 
              end
              
              if isempty(cc21)
                  c21 = eps;
              else
                  c21 = mean(L(cc21));
              end
              
              if isempty(cc22)
                  c22 = eps;
              else
                  c22 = mean(L(cc22));
              end
                          
              
              f_image11=(L-c11).^2.*Heaviside(phi1).*Heaviside(phi2)+f_image11;
              f_image12=(L-c12).^2.*Heaviside(phi1).*(1-Heaviside(phi2))+f_image12;
              f_image21=(L-c21).^2.*(1-Heaviside(phi1)).*Heaviside(phi2)+f_image21;
              f_image22=(L-c22).^2.*(1-Heaviside(phi1)).*(1-Heaviside(phi2))+f_image22;
              end
              
                 
              curvature1 = mu*kappa(phi1);
              curvature1 = curvature1(nb1);
              
              fim1 = 1/layer.*(-f_image11(nb1)+f_image21(nb1)-f_image12(nb1)+f_image22(nb1));
              fim1 = fim1./max(abs(fim1)+eps);

              
              curvature2 = mu*kappa(phi2);
              curvature2 = curvature2(nb2);
             
              fim2 = 1/layer.*(-f_image11(nb2)+f_image12(nb2)-f_image21(nb2)+f_image22(nb2));
              fim2 = fim2./max(abs(fim2)+eps);

             
              force1 = curvature1+fim1;
              force2 = curvature2+fim2;
             
              dt = 1.5;
              
              old(:,:,1) = phi1;
              old(:,:,2) = phi2;

              
              phi1(nb1) = phi1(nb1)+dt.*force1;
              phi2(nb2) = phi2(nb2)+dt.*force2;
              
              new(:,:,1) = phi1;
              new(:,:,2) = phi2;
              
              indicator = checkstop(old,new,dt);

              if indicator 
                 showphi(I, new, n);
                 
                 seg11 = (phi1>0 & phi2>0); 
                 seg12 = (phi1>0 & phi2<0);
                 seg21 = (phi1<0 & phi2>0);
                 seg22 = (phi1<0 & phi2<0);

                 se = strel('disk',1);
                 aa1 = imerode(seg11,se);
                 aa2 = imerode(seg12,se);
                 aa3 = imerode(seg21,se);
                 aa4 = imerode(seg22,se);
                 seg = aa1+2*aa2+3*aa3+4*aa4;
                
                  figure();imagesc(seg);axis image;title('Final OD Segmentation');

                  return
              end
             
              phi1 = reinitialization(phi1, 0.6);
              phi2 = reinitialization(phi2, 0.6);
             
              if(mod(n,20) == 0) 
                 phi(:,:,1) = phi1;
                 phi(:,:,2) = phi2;
                 showphi(I, phi, n);
              end;
          end;
          phi(:,:,1) = phi1;
          phi(:,:,2) = phi2;
          showphi(I, phi, n);          
        seg11 = (phi1>0 & phi2>0); 
        seg12 = (phi1>0 & phi2<0);
        seg21 = (phi1<0 & phi2>0);
        seg22 = (phi1<0 & phi2<0);

        se = strel('disk',1);
        aa1 = imerode(seg11,se);
        aa2 = imerode(seg12,se);
        aa3 = imerode(seg21,se);
        aa4 = imerode(seg22,se);
        seg = aa1+2*aa2+3*aa3+4*aa4;
        
     figure();imagesc(seg);axis image;title('Final OD Segmentation');
end

