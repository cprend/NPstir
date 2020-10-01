function [dc,fx,fy]=upw(p,c,dx,Lx,Ly)
  n=size(c,1);
  l=[n,1:n];
  r=[1:n,1];
  u=-1/dx*diff(p);
  %u=u*10;
  cl=[c(:,n)-Lx,c];
  cr=[c,c(:,1)+Lx];
  ca=cr+cl; %staggered grid
  cd=cr-cl; %tracer defined in the middle of cell but fluxes need to be defined
  %at the edge of the cell
  fx=0.5*(u.*ca-abs(u).*cd); %if you're not going with the flow it'll be 0 so like
  %an if statement of upwind 
  v=1/dx*diff(p,1,2);
  %v=v*10;
  cl=[c(n,:)-Ly;c];
  cr=[c;c(1,:)+Ly];
  ca=cl+cr;
  cd=cr-cl;
  fy = 0.5*(v.*ca-abs(v).*cd);
  dc=-1/dx*(diff(fx,1,2)+diff(fy)); %divergence of flux
