function dc=diffuse(c,k,dx,Lx,Ly)
  n=size(c,1);np1=n+1;np2=n+2;
  ce=[c(:,n)-Lx,c,c(:,1)+Lx];
  ce=[ce(n,:)-Ly;ce;ce(1,:)+Ly];
  dc=-k/dx^2*(-4*c+ce(2:np1,1:n)+ce(2:np1,3:np2)+ce(1:n,2:np1)+ce(3:np2,2:np1));
