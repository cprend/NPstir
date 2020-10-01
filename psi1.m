function p=psi1(xp,yp,kl,a,th)
global amp;
p=0;
for j=1:size(a,1)
  p=p+amp*a(j)*cos(kl(j,1)*xp+kl(j,2)*yp+th(j));
end
