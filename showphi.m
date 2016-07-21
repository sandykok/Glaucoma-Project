function showphi(I, phi, i)

for j = 1:size(phi,3)
    phi_{j} = phi(:,:,j);
end
  imshow(I,'initialmagnification','fit','displayrange',[0 255]);
  hold on;

  if size(phi,3) == 1
      contour(phi_{1}, [0 0], 'r','LineWidth',4);
      contour(phi_{1}, [0 0], 'g','LineWidth',1.3);
  else
      contour(phi_{1}, [0 0], 'r','LineWidth',4);
      contour(phi_{1}, [0 0], 'x','LineWidth',1.3);
      contour(phi_{2}, [0 0], 'g','LineWidth',4);
      contour(phi_{2}, [0 0], 'x','LineWidth',1.3);
  end
  hold off; 
  title([num2str(i) ' Iterations']); 
  drawnow;