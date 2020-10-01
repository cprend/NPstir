%%%Coupled physical-biological model used in Prend et al. (submitted, GRL). The biological model is a 2 component nutrient-phytoplankton system governed by Equations (2.1.1-2.1.2) and the physical model is a 2-D stirring flow given by Equation (2.2.1). The script calls 3 functions, diffuse.m, psi1.m, and upw.m   

global amp
nx=256; %setting grid size
%setting factors for flow 
dt=1/32; %time stepping
amp=0.1; %velocity amplitude 
dth=0.1*sqrt(dt); %phase
afac=sqrt(12*dt); %factor for streamfunction

%setting model domain and grid
L=4*pi; 
dx=L/nx;
x=[0.5:nx]*dx; 
y=[0.5:nx]'*dx;

[xg,yg]=meshgrid(x,y); %grid for tracer
[xp,yp]=meshgrid([0:nx]*dx); %grid for streamfunction

%wavenumbers for streamfunction
kl=[
0.5,0
0,0.5
1,1
1,-1
];
a=rand(size(kl,1),1);

th=2*pi*rand(size(a,1),1); %phase
agamma=1/20; %factor for streamfunction
kappa=1e-4; %if including numerical diffusion 


%setting initial distributions of N and P  
t=0; %start from time zero 
tmax=1000; %how many timesteps to run for

N00=0.8;N01=0.4*N00;%set initial N
mu0=0.08;mu1=0.4*mu0;%set mu
la=0.03;%set entrainment rate lambda
    
k=1/2;
S0=N00-N01*cos(k*xg);%setting initial S0 distribution
mu=mu0-mu1*cos(k*yg);%setting mu distribution
P0=S0-la./mu;%setting initial P distribution
N0=S0-P0;%setting initial N distribution

%initialize the tracer fields
N=N0;
P=P0;
S=S0;

%run the model
while t<=tmax
  if(rem(t,1)==0)
    imagesc(x,y,[N,P]); %plot N and P at each timestep
    axis('xy','equal')
    title(sprintf("t = %3f",t));
    set(gcf, 'color', 'w');
    set(gca, 'color', 'w');
    colorbar();
    drawnow();
  end
  p=psi1(xp,yp,kl,a,th); %define the streamfunction at each time step using psi1.m function
  [dc1,fx1,fy1]=upw(p,N,dx,0,0); %advect N and P by first order upwind advection scheme using upw.m function 
  [dc2,fx2,fy2]=upw(p,P,dx,0,0);
  [dc3,fx3,fy3]=upw(p,S,dx,0,0);
  
  %apply the advection, diffusion, and biological reactions to find the tracer evolution at each time step
  dNdt=dc1+diffuse(N,kappa,dx,0,0)-mu.*N.*P-la*(N-S0);
  dPdt=dc2+diffuse(P,kappa,dx,0,0)+mu.*N.*P-la*P; 
  dSdt=dc3+diffuse(S,kappa,dx,0,0)+la*(S-S0);
  N=N+dt*dNdt;
  P=P+dt*dPdt;
  S=S+dt*dSdt;
  t=t+dt;
  
  %phase shift and amplitude for renovating wave model
  th=th+dth*(rand(size(a))-0.5);
  a=a+afac*(rand(size(a))-0.5)-agamma*a*dt;
  a=a.*(a<=5)+5*(a>5);
  a=a.*(a>=-5)-5*(a<-5);
end