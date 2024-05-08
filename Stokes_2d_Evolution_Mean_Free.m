%This program solves the 2d Stokes using Backward Euler and mean free pressure basis functions.

%Is Working

clear
close all

%Define the constant values for initial time, end time
nu = 1;
t0 = 0;
Tend = 1;
vep = 0.0001; %Singular limit parameter
meas = 1; %Measure of the geomtery \calO

%% Load the mesh (this loads in all the relevant mesh based parameters). 

Stokes_Mesh

%% Compute the appropriate time step.

% Compute vedge, edge as a vector, and area of each element
ve(:,:,1) = msh.POS(msh.TRIANGLES6(:,3),:) - msh.POS(msh.TRIANGLES6(:,2),:);
ve(:,:,2) = msh.POS(msh.TRIANGLES6(:,1),:) - msh.POS(msh.TRIANGLES6(:,3),:);
ve(:,:,3) = msh.POS(msh.TRIANGLES6(:,2),:) - msh.POS(msh.TRIANGLES6(:,1),:);
AREA = 0.5*abs(-ve(:,1,3).*ve(:,2,2) + ve(:,2,3).*ve(:,1,2));

dtf = 0.1*(min(AREA))/nu; %CFL condition is ~0.5*area/nu
% dtf = 10*dtf;

time_steps = round((Tend-t0)/dtf);

%Ensures that t0 + dtf*time_steps = Tend;
dtf = (Tend - t0)/time_steps;

clear ve AREA


%% Define the standard fluid element
phi1 = @(x,y) (1-x-y)*(1-2*x-2*y);
phi2 = @(x,y) y*(2*y-1);
phi3 = @(x,y) x*(2*x-1);
phi4 = @(x,y) 4*(1-x-y)*y;
phi5 = @(x,y) 4*x*y;
phi6 = @(x,y) 4*(1-x-y)*x;

phi1x = @(x,y) -3+4*x+4*y;
phi2x = @(x,y) 0;
phi3x = @(x,y) 4*x-1;
phi4x = @(x,y) -4*y;
phi5x = @(x,y) 4*y;
phi6x = @(x,y) 4-8*x-4*y;

phi1y = @(x,y) -3+4*x+4*y;
phi2y = @(x,y) 4*y-1;
phi3y = @(x,y) 0;
phi4y = @(x,y) 4-4*x-8*y;
phi5y = @(x,y) 4*x;
phi6y = @(x,y) -4*x;

%Integration points for STRANG5 (see
%https://people.sc.fsu.edu/~jburkardt/datasets/quadrature_rules_tri/quadrature_rules_tri.html)
FLUID_INTPTS = [0.816847572980459,  0.091576213509771;
  0.091576213509771,  0.816847572980459;
  0.091576213509771,  0.091576213509771;
  0.108103018168070,  0.445948490915965;
  0.445948490915965,  0.108103018168070;
  0.445948490915965,  0.445948490915965];

%Integration weights for STRANG5
 FLUID_WTS = [0.109951743655322, 0.109951743655322, 0.109951743655322,0.223381589678011,0.223381589678011,0.223381589678011];

%Evaluates phi at the integration points.
nbFluidintpts = length(FLUID_INTPTS(:,1));

phi = zeros(6, nbFluidintpts);
phix = zeros(6, nbFluidintpts);
phiy = zeros(6, nbFluidintpts);

for i = 1:1:nbFluidintpts
    phi(1,i) = phi1(FLUID_INTPTS(i,1),FLUID_INTPTS(i,2));
    phi(2,i) = phi2(FLUID_INTPTS(i,1),FLUID_INTPTS(i,2));
    phi(3,i) = phi3(FLUID_INTPTS(i,1),FLUID_INTPTS(i,2));
    phi(4,i) = phi4(FLUID_INTPTS(i,1),FLUID_INTPTS(i,2));
    phi(5,i) = phi5(FLUID_INTPTS(i,1),FLUID_INTPTS(i,2));
    phi(6,i) = phi6(FLUID_INTPTS(i,1),FLUID_INTPTS(i,2));
    
    phix(1,i) = phi1x(FLUID_INTPTS(i,1),FLUID_INTPTS(i,2));
    phix(2,i) = phi2x(FLUID_INTPTS(i,1),FLUID_INTPTS(i,2));
    phix(3,i) = phi3x(FLUID_INTPTS(i,1),FLUID_INTPTS(i,2));
    phix(4,i) = phi4x(FLUID_INTPTS(i,1),FLUID_INTPTS(i,2));
    phix(5,i) = phi5x(FLUID_INTPTS(i,1),FLUID_INTPTS(i,2));
    phix(6,i) = phi6x(FLUID_INTPTS(i,1),FLUID_INTPTS(i,2));
    
    phiy(1,i) = phi1y(FLUID_INTPTS(i,1),FLUID_INTPTS(i,2));
    phiy(2,i) = phi2y(FLUID_INTPTS(i,1),FLUID_INTPTS(i,2));
    phiy(3,i) = phi3y(FLUID_INTPTS(i,1),FLUID_INTPTS(i,2));
    phiy(4,i) = phi4y(FLUID_INTPTS(i,1),FLUID_INTPTS(i,2));
    phiy(5,i) = phi5y(FLUID_INTPTS(i,1),FLUID_INTPTS(i,2));
    phiy(6,i) = phi6y(FLUID_INTPTS(i,1),FLUID_INTPTS(i,2));
end

clear phi1 phi2 phi3 phi4 phi5 phi6
clear phi1x phi2x phi3x phi4x phi5x phi6x
clear phi1y phi2y phi3y phi4y phi5y phi6y

%% Define the standard pressure element
psi1 = @(x,y) 1-x-y;
psi2 = @(x,y) y;
psi3 = @(x,y) x;

psi = zeros(3,nbFluidintpts);

for i = 1:1:nbFluidintpts
    psi(1,i) = psi1(FLUID_INTPTS(i,1),FLUID_INTPTS(i,2));
    psi(2,i) = psi2(FLUID_INTPTS(i,1),FLUID_INTPTS(i,2));
    psi(3,i) = psi3(FLUID_INTPTS(i,1),FLUID_INTPTS(i,2));
end

clear psi1 psi2 psi3

%% Computes the mass, stiffness, and B matrices for fluid and pressure.
kf = sparse(Mf,Mf);
Massf = sparse(Mf,Mf);
Bx = sparse(Mf,Mp);
By = sparse(Mf,Mp);
Massp = sparse(Mp,Mp);
Ep = zeros(1,Mp); %Projections onto mean free space
Efx = zeros(1,Mf); %Integrals of phix for mean free div matrix
Efy = zeros(1,Mf); %Integrals of phiy for mean free div matrix

for ell = 1:1:Lf
    x_1 = msh.POS(msh.TRIANGLES6(ell,1),1);
    y_1 = msh.POS(msh.TRIANGLES6(ell,1),2);
    
    x_2 = msh.POS(msh.TRIANGLES6(ell,2),1);
    y_2 = msh.POS(msh.TRIANGLES6(ell,2),2);
    
    x_3 = msh.POS(msh.TRIANGLES6(ell,3),1);
    y_3 = msh.POS(msh.TRIANGLES6(ell,3),2);
    
    D = (x_3 - x_1)*(y_2 - y_1) - (x_2 - x_1)*(y_3 - y_1);
    
    sgn = sign(D);

    E_1 = (x_2 - x_1)^2 + (y_2 - y_1)^2;
    E_2 = (x_3 - x_1)*(x_2 - x_1) + (y_3 - y_1)*(y_2 - y_1);
    E_3 = (x_3 - x_1)^2 + (y_3 - y_1)^2;

    xtildex =  (y_2-y_1);
    xtildey = -(x_2-x_1);
    ytildex = -(y_3-y_1);
    ytildey =  (x_3-x_1);
    
    for m = 1:3
        j = msh.TH3(ell,m);
        Ep(j) = Ep(j) + abs(D)/2 * sum(FLUID_WTS.*psi(m,:));
        for n = 1:3
            i = msh.TH3(ell,n);
            Massp(i,j) = Massp(i,j) + 0.5*abs(D)*sum(FLUID_WTS.*psi(n,:).*psi(m,:));
        end
    end

    for r = 1:1:Tf
        i = msh.TRIANGLES6(ell,r); 
        
        Efx(i) = Efx(i) + sgn/2 * sum(FLUID_WTS.*phix(r,:));
        Efy(i) = Efy(i) + sgn/2 * sum(FLUID_WTS.*phiy(r,:));
        
        for s = 1:1:Tf
            k1 = E_1*sum(FLUID_WTS.*phix(r,:).*phix(s,:));
            k2 = E_2*(sum(FLUID_WTS.*phix(r,:).*phiy(s,:)) + sum(FLUID_WTS.*phiy(r,:).*phix(s,:)));
            k3 = E_3*sum(FLUID_WTS.*phiy(r,:).*phiy(s,:));
            
            k_tilde = 0.5*abs(D)^(-1)*(k1-k2+k3);
            massf_tilde = 0.5*abs(D)*sum(FLUID_WTS.*phi(r,:).*phi(s,:));
           
            j = msh.TRIANGLES6(ell,s); 
            
            kf(i,j) = kf(i,j) + k_tilde;
            Massf(i,j) = Massf(i,j) + massf_tilde;
        end
         for m = 1:1:3
             bxtilde = -0.5*sgn*sum(FLUID_WTS.*psi(m,:).*(xtildex*phix(r,:)+ytildex*phiy(r,:)));
             bytilde = -0.5*sgn*sum(FLUID_WTS.*psi(m,:).*(xtildey*phix(r,:)+ytildey*phiy(r,:)));
             
             j = msh.TH3(ell,m);
            
             Bx(i,j) = Bx(i,j) + bxtilde;
             By(i,j) = By(i,j) + bytilde;
         end
    end
    

end

%This, along with the \vep Mass_p, is the culprit of the memory errors.
Bx = sparse(Bx - 1/meas*transpose(Efx)*Ep);
By = sparse(By - 1/meas*transpose(Efy)*Ep);

clear k1 k2 k3 k_tilde massf_tilde bxtilde bytilde D E_1 E_2 E_3 sgn x_1 x_2 x_3 y_1 y_2 y_3 xtildex xtildey ytildex ytildey



%% Define the initial data

x = msh.POS(:,1); %x coordinates of fluid nodes in CalO
y = msh.POS(:,2); %y coordinates of fluid nodes in CalO

%Non-zero boundary data, nonzero variable forcing

%True solutions: v1 = exp(-t)sin^2(pi x)sin(2 pi y)
% v2 = exp(-t)sin(2 pi x)sin^2(pi y)
% p = (1 + exp(-t))(x - y)


%Initial fluid data
vn1 = (sin(pi*x)).^2.*sin(2*pi*y); %Initial x component of velocity
vn2 = -sin(2*pi*x).*(sin(pi*y)).^2; %Initial y component of velocity

%Saves the initial datat to be plotted later.
v01 = vn1;
v02 = vn2;


xp = msh.POS(Jp,1);
yp = msh.POS(Jp,2);

%Per Ern & Geurmond Theorem 6.43, this initial pressure data is arbitrary
%unless you want a more accurate in time method.
%p0 = 2*(xp - yp); %"True initial data"

%p0 = 7*xp; %A random guess
p0 = rand(Mp,1); %An even more random guess.

pn = p0;

clear x y xp yp


%% Builds the fluid and pressure forcing vectors. Forcing is taken to be of the form f1(t)*Ff + f2(t)*Fp, where Ff and Fp do not depend on time.


f1pvals = ones(nbFluidintpts,1);
f2pvals = -f1pvals;


F1f = zeros(Mf,1);
F2f = zeros(Mf,1);

F1p = zeros(Mf,1);
F2p = zeros(Mf,1);

Fdiv = zeros(Mp,1); %The forcing associated to divergence. Will be zero in final program.

for ell = 1:1:Lf
    x_1 = msh.POS(msh.TRIANGLES6(ell,1),1);
    y_1 = msh.POS(msh.TRIANGLES6(ell,1),2);
    
    x_2 = msh.POS(msh.TRIANGLES6(ell,2),1);
    y_2 = msh.POS(msh.TRIANGLES6(ell,2),2);
    
    x_3 = msh.POS(msh.TRIANGLES6(ell,3),1);
    y_3 = msh.POS(msh.TRIANGLES6(ell,3),2);
    
    D = (x_3 - x_1)*(y_2 - y_1) - (x_2 - x_1)*(y_3 - y_1);
    
   
    glbpts = [x_3-x_1, x_2-x_1; y_3-y_1, y_2-y_1]*transpose(FLUID_INTPTS) + [x_1;y_1];
    glbpts = transpose(glbpts);
      

    %Forcing associated to fluid
     x = glbpts(:,1);
     y = glbpts(:,2);
     f1fvals = -2*pi^2*cos(pi*x).^2.*sin(2*pi*y) + 6*pi^2*sin(pi*x).^2.*sin(2*pi*y) - sin(pi*x).^2.*sin(2*pi*y);
     f2fvals = 2*pi^2*cos(pi*y).^2.*sin(2*pi*x) - 6*pi^2*sin(2*pi*x).*sin(pi*y).^2 + sin(2*pi*x).*sin(pi*y).^2;
    
     %Divergence forcing
     fdiv_vals = x - y;

 for r = 1:1:Tf
        i = msh.TRIANGLES6(ell,r); 
          
        F1f_tilde = abs(D)/2 * sum(FLUID_WTS.*phi(r,:)*f1fvals);
        F2f_tilde = abs(D)/2 * sum(FLUID_WTS.*phi(r,:)*f2fvals);
        F1f(i) = F1f(i) + F1f_tilde;
        F2f(i) = F2f(i) + F2f_tilde;
        
        F1p_tilde = abs(D)/2 * sum(FLUID_WTS.*phi(r,:)*f1pvals);
         F2p_tilde = abs(D)/2 * sum(FLUID_WTS.*phi(r,:)*f2pvals);
        F1p(i) = F1p(i) + F1p_tilde;
         F2p(i) = F2p(i) + F2p_tilde;
         
      
 end
 
 for r = 1:3
     i = msh.TH3(ell,r);
     
      Fdiv_tilde = abs(D)/2 * sum(FLUID_WTS.*psi(r,:)*fdiv_vals);
      Fdiv(i) = Fdiv(i) + Fdiv_tilde;
 end
end


clear Fdiv_tilde fdiv_vals F1f_tilde F2f_tilde f1fvals f2fvals F1p_tilde F2p_tilde f1pvals f2pvals D x_1 x_2 x_3 y_1 y_2 y_3 x y


%% Arrange the LHSf as a Block Matrix

%Fluid LHS
LHSf = sparse(2*Mf + Mp,2*Mf + Mp);

LHSf(1:Mf,1:Mf) = Massf + dtf*nu*kf; %11 block
LHSf(Mf+1:2*Mf,Mf+1:2*Mf) = Massf + dtf*nu*kf; %22 block
LHSf(1:Mf,2*Mf+1:2*Mf+Mp) = dtf*Bx; %13 block
LHSf(Mf+1:2*Mf,2*Mf+1:2*Mf+Mp) = dtf*By; %23 block
LHSf(2*Mf+1:2*Mf+Mp,1:Mf) = -dtf*transpose(Bx); %31 block
LHSf(2*Mf+1:2*Mf+Mp,Mf+1:2*Mf) = -dtf*transpose(By); %32 block
LHSf(2*Mf+1:2*Mf+Mp,2*Mf+1:2*Mf+Mp) = vep*Massp;


%% Loops through fluid nodes to ensure Dirichlet boundary conditions are satisfied

   LHSf(J0f,:) = 0;
   LHSf(J0f + Mf,:) = 0;
    
   LHSf(1:Mf,J0f) = 0;
   LHSf(Mf+1:2*Mf,J0f + Mf) = 0;


   LHSf(J0f,J0f) = eye(M0f);
   LHSf(J0f + Mf,J0f + Mf) = eye(M0f);

%% Time stepping fluid
% A = Massf + nu*dtf*kf;

RHSf = zeros(2*Mf + Mp,1);


    
    dLHSf = decomposition(LHSf);
    
%    clear LHSf
    
 time = t0 + dtf;
    
    for tk = 2:1:time_steps
    RHSf(1:Mf) = dtf*exp(-time)*F1f + dtf*(1+exp(-time))*F1p + Massf*vn1;
    RHSf(Mf+1:2*Mf) = dtf*exp(-time)*F2f + dtf*(1+exp(-time))*F2p + Massf*vn2; 
    RHSf(2*Mf+1:2*Mf + Mp) = vep*dtf*(-exp(-time))*Fdiv + vep*Massp*pn;

   %Makes the RHS vector satisfy the Dirichlet data
   %This section handles interior nodes.
%    RHSf(Jf + Mf) = RHSf(Jf + Mf) - A(Jf,J0f)*gamma2t(J0f);
   RHSf(J0f) = 0;
%    RHSf(J0f + Mf) = gamma2t(J0f);
   RHSf(J0f + Mf) = 0;

   
   %Solve for the coefficients of fluids and pressure at time
   %step N+1
   coefficients_f = full(dLHSf\RHSf);
   vn1 = coefficients_f(1:Mf);
   vn2 = coefficients_f(Mf+1:2*Mf);
   pn = coefficients_f(2*Mf+1:2*Mf+Mp);
   time = time + dtf;
    end

clear A coefficients_f
%% Plots the FEM and True Solutions

plot_soln = 1; %Plot the solution? 1 = yes, 0 = no.
plot_mesh = 0; %Plot the mesh? 1 = yes, 0 = no.
 
if plot_mesh == 1
    meshplotter2D
end

if plot_soln == 1
   
    x = msh.POS(:,1);
    y = msh.POS(:,2);
    xp = msh.POS(Jp,1);
    yp = msh.POS(Jp,2);
    
    %Initial Data
%    v01 = 6*(y.^2 - y).*(x-1).^3.*x.^3; %Initial x component of velocity
%    v02 = -3*(2*y.^3 - 3*y.^2).*x.^2.*(x-1).^2.*(2*x-1); %Initial y component of velocity
%    p0 = 2*(xp - yp);
    
   v1T = exp(-Tend)*(sin(pi*x)).^2.*sin(2*pi*y); %Final x component of velocity
   v2T = -exp(-Tend)*sin(2*pi*x).*(sin(pi*y)).^2; %Final y component of velocity
   pT = (1+exp(-Tend))*(xp - yp);
   

       figure(4)
        trisurf(msh.TRIANGLES6(:,1:3), msh.POS(:,1),msh.POS(:,2), vn1);
        shading interp;
        axis('square');
        xlabel('x','fontsize',14)
        ylabel('y','fontsize',14)
%         view(2);
         title('FEM X Velocity');

          figure(5)
        trisurf(msh.TRIANGLES6(:,1:3), msh.POS(:,1),msh.POS(:,2), v1T);
        shading interp;
        axis('square');
        xlabel('x','fontsize',14)
        ylabel('y','fontsize',14)
%         view(2);
         title('X Velocity');   
         
         figure(6)
        trisurf(msh.TRIANGLES6(:,1:3), msh.POS(:,1),msh.POS(:,2), abs(v1T-vn1));
        shading interp;
        axis('square');
        xlabel('x','fontsize',14)
        ylabel('y','fontsize',14)
        zlim([0, max(abs(v1T))])
        caxis([0, max(abs(v1T))])
%         view(2);
         title('X Velocity minus FEM X velocity');  
         

     figure(7)
        trisurf(msh.TRIANGLES6(:,1:3), msh.POS(:,1),msh.POS(:,2), vn2);
        shading interp;
        axis('square');
        % axis([0,1,0,1]);
        xlabel('x','fontsize',14)
        ylabel('y','fontsize',14)
%         view(2);
         title('FEM Y Velocity');

              figure(8)
        trisurf(msh.TRIANGLES6(:,1:3), msh.POS(:,1),msh.POS(:,2), v2T);
        shading interp;
        axis('square');
        % axis([0,1,0,1]);
        xlabel('x','fontsize',14)
        ylabel('y','fontsize',14)
%         view(2);
         title('Y Velocity');
  
        figure(9)
        trisurf(msh.TRIANGLES6(:,1:3), msh.POS(:,1),msh.POS(:,2), abs(v2T-vn2));
        shading interp;
        axis('square');
        % axis([0,1,0,1]);
        xlabel('x','fontsize',14)
        ylabel('y','fontsize',14)
        zlim([0, max(abs(v2T))])
        caxis([0, max(abs(v2T))])
%         view(2);
         title('Y Velocity minus FEM Y velocity');

         
   figure(10)
        trisurf(msh.TH3(:,1:3), msh.POS(Jp,1),msh.POS(Jp,2), pn);
        shading interp;
        axis('square');
        % axis([0,1,0,1]);
        xlabel('x','fontsize',14)
        ylabel('y','fontsize',14)
%         view(2);
         title('FEM Pressure');

        figure(11)
        trisurf(msh.TH3(:,1:3), msh.POS(Jp,1),msh.POS(Jp,2), pT);
        shading interp;
        axis('square');
        % axis([0,1,0,1]);
        xlabel('x','fontsize',14)
        ylabel('y','fontsize',14)
%       view(2);
        title('True Pressure');
        
        figure(12)
        trisurf(msh.TH3(:,1:3), msh.POS(Jp,1),msh.POS(Jp,2), abs(pT-pn));
        shading interp;
        axis('square');
        % axis([0,1,0,1]);
        xlabel('x','fontsize',14)
        ylabel('y','fontsize',14)
        zlim([0, max(abs(pT))])
        caxis([0, max(abs(pT))])
%       view(2);
        title('True Pressure minus FEM Pressure');
         
end         

clear plot_soln x y xp yp plot_mesh
   