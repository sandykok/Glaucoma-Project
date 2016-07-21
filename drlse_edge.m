function phi = drlse_edge(phi_0, g, lambda,mu, alfa, epsilon, timestep, iter, potentialFunction)
phi=phi_0;
[vx, vy]=gradient(g);
for k=1:iter
    phi=NeumannBoundCond(phi);
    [phi_x,phi_y]=gradient(phi);
    s=sqrt(phi_x.^2 + phi_y.^2);
    smallNumber=1e-10;  
    Nx=phi_x./(s+smallNumber); 
    Ny=phi_y./(s+smallNumber);
    curvature=div(Nx,Ny);
    if strcmp(potentialFunction,'single-well')
        distRegTerm = 4*del2(phi)-curvature;  
    elseif strcmp(potentialFunction,'double-well');
        distRegTerm=distReg_p2(phi); 
    else
        disp('Error: Wrong choice of potential function. Please input the string "single-well" or "double-well" in the drlse_edge function.');
    end        
    diracPhi=Dirac(phi,epsilon);
    areaTerm=diracPhi.*g; 
    edgeTerm=diracPhi.*(vx.*Nx+vy.*Ny) + diracPhi.*g.*curvature;
    phi=phi + timestep*(mu*distRegTerm + lambda*edgeTerm + alfa*areaTerm);
end


function f = distReg_p2(phi)
[phi_x,phi_y]=gradient(phi);
s=sqrt(phi_x.^2 + phi_y.^2);
a=(s>=0) & (s<=1);
b=(s>1);
ps=a.*sin(2*pi*s)/(2*pi)+b.*(s-1);  
dps=((ps~=0).*ps+(ps==0))./((s~=0).*s+(s==0));
f = div(dps.*phi_x - phi_x, dps.*phi_y - phi_y) + 4*del2(phi);  

function f = div(nx,ny)
[nxx,junk]=gradient(nx);  
[junk,nyy]=gradient(ny);
f=nxx+nyy;

function f = Dirac(x, sigma)
f=(1/2/sigma)*(1+cos(pi*x/sigma));
b = (x<=sigma) & (x>=-sigma);
f = f.*b;

function g = NeumannBoundCond(f)
[nrow,ncol] = size(f);
g = f;
g([1 nrow],[1 ncol]) = g([3 nrow-2],[3 ncol-2]);  
g([1 nrow],2:end-1) = g([3 nrow-2],2:end-1);          
g(2:end-1,[1 ncol]) = g(2:end-1,[3 ncol-2]);  