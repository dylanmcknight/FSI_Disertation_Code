%This program solves the stationary spectral 2d-1d Fluid structure interaction problem

%Is working
%Requires:Thesis_FSI_Mesh_Quintic.m, Spectral_Norm.m

clear
close all

trials = 100;
Solnnorm = zeros(trials,1);
Fnorm = zeros(trials,1);

for blah = 1:1:trials

%Define the constant values for beta
nu = 1;
meas = 1;
beta = blah*1i; %\beta \Phi - \scrA \Phi = \Phi*

penalty = 10^4; %parameter to enforce boundary conditions for fluid, may have to tune for othe problems
tol = 10^(-14); %parameter to clean up mean enformxing vector
%Both of these parameters are surely mesh dependent, hence the need to
%tune. The above values are good for the given test problems and mesh
%files.


%% Load the mesh (this loads in all the relevant mesh based parameters). 

Thesis_FSI_Mesh_Quintic


%% Define the standard plate element.

theta1 = @(x) -17/128 + 21/128*x + 81/64*x^2 -101/64*x^3 -81/128*x^4 + 117/128*x^5; %w(-1)
theta2 = @(x) -17/128 - 21/128*x + 81/64*x^2 +101/64*x^3 -81/128*x^4 - 117/128*x^5; %w(1)
theta3 = @(x) 81/128 -243/128*x -81/64*x^2 +243/64*x^3 +81/128*x^4 -243/128*x^5; %w(-1/3)
theta4 = @(x) 81/128 +243/128*x -81/64*x^2 -243/64*x^3 +81/128*x^4 +243/128*x^5; %w(1/3)
theta5 = @(x) -1/32 + 1/32*x + 5/16*x^2 -5/16*x^3 -9/32*x^4 +9/32*x^5; %w'(-1)
theta6 = @(x) 1/32 + 1/32*x - 5/16*x^2 -5/16*x^3 +9/32*x^4 +9/32*x^5; %w'(1)                    

theta1xx = @(x) 585/32*x^3 - 243/32*x^2 - 303/32*x + 81/32; %w(-1)
theta2xx = @(x) -585/32*x^3 - 243/32*x^2 + 303/32*x + 81/32; %w(1)
theta3xx = @(x) -1215/32*x^3 + 243/32*x^2 + 729/32*x - 81/32; %w(-1/3)
theta4xx = @(x) 1215/32*x^3 + 243/32*x^2 - 729/32*x - 81/32; %w(1/3)
theta5xx = @(x) 45/8*x^3 - 27/8*x^2 - 15/8*x + 5/8; %w'(-1)
theta6xx = @(x) +45/8*x^3 + 27/8*x^2 - 15/8*x - 5/8; %w'(1)

%Using local node ordering [1)--(3)--(4)--(2]

%Gaussian Quadrature on [-1, 1], see Abramowitz & Stegun
% PLATE_INTPTS = [-0.906179845938664,-0.538469310105683,0,0.538469310105683,0.906179845938664];
% PLATE_WTS = [0.236926885056189,0.478628670499366,0.568888888888889,0.478628670499366,0.236926885056189];

PLATE_INTPTS = [-0.932469514203152, -0.661209386466265, -0.238619186083197, 0.238619186083197, 0.661209386466265, 0.932469514203152];
PLATE_WTS = [0.171324492379170, 0.360761573048139, 0.467913934572691, 0.467913934572691, 0.360761573048139, 0.171324492379170];


nbPlateintpts = length(PLATE_INTPTS);


theta = zeros(6, nbPlateintpts);
for i = 1:1:nbPlateintpts
    theta(:,i) = [theta1(PLATE_INTPTS(i)) ; theta2(PLATE_INTPTS(i)); theta3(PLATE_INTPTS(i)); theta4(PLATE_INTPTS(i)); theta5(PLATE_INTPTS(i)); theta6(PLATE_INTPTS(i))];
end

thetaxx = zeros(6, nbPlateintpts);
for i = 1:1:nbPlateintpts
    thetaxx(:,i) = [theta1xx(PLATE_INTPTS(i)) ; theta2xx(PLATE_INTPTS(i)); theta3xx(PLATE_INTPTS(i)); theta4xx(PLATE_INTPTS(i)); theta5xx(PLATE_INTPTS(i)); theta6xx(PLATE_INTPTS(i))];
end

%thetamid evaluates the structure basis functions at the midpoint, as we
%need this value when computing the boundary data matching.
thetamid = [theta1(0); theta2(0); theta3(0); theta4(0); theta5(0); theta6(0)];

clear theta1 theta2 theta3 theta4 theta5 theta6
clear theta1xx theta2xx theta3xx theta4xx theta5xx theta6xx

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

%% Computes the mass, stiffness, B matrices, and mean enforcing vector for fluid and pressure.
kf = sparse(Mf,Mf);
Massf = sparse(Mf,Mf);
Bx = sparse(Mf,Mp);
By = sparse(Mf,Mp);
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

for i = 1:Mf
 if abs(Efx(i)) < tol
     Efx(i) = 0;
 end
end


for i = 1:Mf
 if abs(Efy(i)) < tol
     Efy(i) = 0;
 end
end

for i = 1:Mp
    if abs(Ep(i)) < tol
        Ep(i) = 0;
    end
end


% Bx = Bx - 1/meas*transpose(Efx)*Ep;
% By = By - 1/meas*transpose(Efy)*Ep;


clear k1 k2 k3 k_tilde massf_tilde bxtilde bytilde D E_1 E_2 E_3 sgn x_1 x_2 x_3 y_1 y_2 y_3 xtildex xtildey ytildex ytildey

%% Computes the plate mass and S matrices
S = sparse(dofs,dofs);
Mass_s = sparse(dofs,dofs);
Es = zeros(1,dofs);

for ell = 1:1:Ls   
    %Handles the Hermite Endpoints
     for r = 1:1:2
        i = msh.LOCSTRUCT(ell,r);
         
         Es(i) = Es(i) + (0.5*h)*sum(PLATE_WTS.*theta(r,:));
         Es(i + Ms) = Es(i + Ms) + (0.5*h)^2*sum(PLATE_WTS.*theta(r + 4,:));
         
        for s = 1:1:2
            
            j = msh.LOCSTRUCT(ell,s);
            
            S(i,j) = S(i,j) + 1/(0.5*h)^3*sum(PLATE_WTS.*thetaxx(r,:).*thetaxx(s,:)); %w with w
            S(i + Ms, j) = S(i + Ms, j) + 1/(0.5*h)^2*sum(PLATE_WTS.*thetaxx(r + 4,:).*thetaxx(s,:)); %w' with w
            S(i,j + Ms) = S(i,j + Ms) + 1/(0.5*h)^2*sum(PLATE_WTS.*thetaxx(r,:).*thetaxx(s + 4,:)); %w with w'
            S(i + Ms, j + Ms) = S(i + Ms, j + Ms) + 1/(0.5*h)*sum(PLATE_WTS.*thetaxx(r + 4,:).*thetaxx(s + 4,:)); %w' with w'
            
            Mass_s(i,j) = Mass_s(i,j) + (0.5*h)*sum(PLATE_WTS.*theta(r,:).*theta(s,:)); %w with w
            Mass_s(i + Ms, j) = Mass_s(i + Ms, j) + (0.5*h)^2*sum(PLATE_WTS.*theta(r + 4,:).*theta(s,:)); %w' with w
            Mass_s(i,j + Ms) = Mass_s(i,j + Ms) + (0.5*h)^2*sum(PLATE_WTS.*theta(r,:).*theta(s + 4,:)); %w with w'
            Mass_s(i + Ms, j + Ms) = Mass_s(i + Ms, j + Ms) + (0.5*h)^3*sum(PLATE_WTS.*theta(r + 4,:).*theta(s + 4,:)); %w' with w'
        end
     end
    %Handles the interaction between Hermite end points and Lagrangian
    %Middle points (r = 3,4, s = 1,2)
    for r = 3:1:4
    i = msh.LOCSTRUCT(ell,r);
    
    Es(i) = Es(i) + (0.5*h)*sum(PLATE_WTS.*theta(r,:));
    
    for s = 1:1:2
       j = msh.LOCSTRUCT(ell,s);
       
       S(i,j) = S(i,j) + 1/(0.5*h)^3*sum(PLATE_WTS.*thetaxx(r,:).*thetaxx(s,:)); %w with w
       S(i,j + Ms) = S(i,j + Ms) + 1/(0.5*h)^2*sum(PLATE_WTS.*thetaxx(r,:).*thetaxx(s + 4,:)); %w with w'
            
       Mass_s(i,j) = Mass_s(i,j) + (0.5*h)*sum(PLATE_WTS.*theta(r,:).*theta(s,:)); %w with w
       Mass_s(i,j + Ms) = Mass_s(i,j + Ms) + (0.5*h)^2*sum(PLATE_WTS.*theta(r,:).*theta(s + 4,:)); %w with w'
    end
    end
    
    %Handles interaction Lagrangian/Lagrangian terms (r = s = 3,4)
    for r = 3:1:4
     i = msh.LOCSTRUCT(ell,r);
        
        for s = 3:1:4
        j = msh.LOCSTRUCT(ell,s);
        S(i,j) = S(i,j) + 1/(0.5*h)^3*sum(PLATE_WTS.*thetaxx(r,:).*thetaxx(s,:));
        Mass_s(i,j) = Mass_s(i,j) + (0.5*h)*sum(PLATE_WTS.*theta(r,:).*theta(s,:));
        end
    end
     
    %Handles the interaction between Hermite end points and Lagrangian
    %Middle points (r = 1,2, s = 3,4)
    for s = 3:1:4
    j = msh.LOCSTRUCT(ell,s);
    for r = 1:1:2
       i = msh.LOCSTRUCT(ell,r);
       
       S(i,j) = S(i,j) + 1/(0.5*h)^3*sum(PLATE_WTS.*thetaxx(r,:).*thetaxx(s,:)); %w with w
       S(i + Ms,j) = S(i + Ms,j) + 1/(0.5*h)^2*sum(PLATE_WTS.*thetaxx(r + 4,:).*thetaxx(s,:)); %w' with w
            
       Mass_s(i,j) = Mass_s(i,j) + (0.5*h)*sum(PLATE_WTS.*theta(r,:).*theta(s,:)); %w with w
       Mass_s(i + Ms,j) = Mass_s(i + Ms,j) + (0.5*h)^2*sum(PLATE_WTS.*theta(r + 4,:).*theta(s,:)); %w' with w
    end
    end
end


%% Builds the fluid and pressure forcing vectors. Forcing is taken to be of the form f1(t)*Ff + f2(t)*Fp, where Ff and Fp do not depend on time.


F1f = zeros(Mf,1);
F2f = zeros(Mf,1);

F1p = zeros(Mf,1);
F2p = zeros(Mf,1);

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
      

    %v2 Nonzero Dirichlet B.C.
     x = glbpts(:,1);
     y = glbpts(:,2);
     f1fvals = 6*beta*1*(y.^2 - y).*(x-1).^3.*x.^3 + nu*(-12*(x - 1).^3.*x.^3 - 36*(y.^2 - y).*(x - 1).^3.*x - 108*(y.^2 - y).*(x - 1).^2.*x.^2 - 36*(y.^2 - y).*(x - 1).*x.^3);
     f2fvals = -3*beta*1*(2*y.^3 - 3*y.^2).*x.^2.*(x-1).^2.*(2*x-1) + nu*(18*(x - 1).^3.*x.^2.*(2*y - 1) + 18*(x - 1).^2.*x.^3.*(2*y - 1) + 6*(2*y.^3 - 3*y.^2).*(x - 1).^3 + 54*(2*y.^3 - 3*y.^2).*(x - 1).^2.*x + 54*(2*y.^3 - 3*y.^2).*(x - 1).*x.^2 + 6*(2*y.^3 - 3*y.^2).*x.^3);
     
     f1pvals = 1;
     f2pvals = -1;

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
end


clear F1f_tilde F2f_tilde f1fvals f2fvals F1p_tilde F2p_tilde f1pvals f2pvals D x_1 x_2 x_3 y_1 y_2 y_3 x y

 %% Constructs the forcing term in w_{tt} + \Aw = F
F1s = zeros(dofs, 1); %F(t,x) can be factored as exp(-t)*f(x) for our test problems. Thus, we will compute F0 = (f(x),theta(x)) and use the explicit time evolution in the time stepping.
F2s = zeros(dofs, 1);
Fc = zeros(dofs,1); %The component of forcing that arises from BB pressure.

for ell = 1:1:Ls
      x1 = msh.LOCPOS(msh.LOCSTRUCT(ell,1),1);
      glbpts = (x1+h/2)*ones(1,nbPlateintpts) + h/2*PLATE_INTPTS; %transformation from [-1,1] to [x1,x2] is x = h/2(tildex) + (x1+h/2)
     
      f1svals = 3*beta*glbpts.^2.*(glbpts-1).^2.*(2*glbpts-1) - 720*glbpts + 360; 
      f2svals = -3*(beta + 1)*glbpts.^2.*(glbpts-1).^2.*(2*glbpts-1);

      fcval = -1; %BB pressure is a constant. It is negative of this number since it is -c in mixed formulation
%Handle the Hermite endpoints
     for r = 1:1:2
       F1s_tildetheta =  (0.5*h) * sum(PLATE_WTS.*f1svals.*theta(r,:));
       F1s_tildethetax =  (0.5*h)^2 * sum(PLATE_WTS.*f1svals.*theta(r+4,:));

       F2s_tildetheta =  (0.5*h) * sum(PLATE_WTS.*f2svals.*theta(r,:));
       F2s_tildethetax =  (0.5*h)^2 * sum(PLATE_WTS.*f2svals.*theta(r+4,:));

       i = msh.LOCSTRUCT(ell,r);
       F1s(i) = F1s(i) + F1s_tildetheta;
       F1s(i+Ms) = F1s(i+Ms) + F1s_tildethetax;
       
       F2s(i) = F2s(i) + F2s_tildetheta;
       F2s(i+Ms) = F2s(i+Ms) + F2s_tildethetax;

       Fc(i) = Fc(i) + (0.5*h)*fcval * sum(PLATE_WTS.*theta(r,:));
       Fc(i + Ms) = Fc(i + Ms) + (0.5*h)^2 * fcval * sum(PLATE_WTS.*theta(r + 4,:));
     end
    %Handle the Lagrangian middle points
    i = msh.LOCSTRUCT(ell,3);
    F1s_tildetheta =  (0.5*h) * sum(PLATE_WTS.*f1svals.*theta(3,:));
    F1s(i) = F1s(i) + F1s_tildetheta;

    F2s_tildetheta =  (0.5*h) * sum(PLATE_WTS.*f2svals.*theta(3,:));
    F2s(i) = F2s(i) + F2s_tildetheta;
    
    Fc(i) = Fc(i) + (0.5*h) * fcval * sum(PLATE_WTS.*theta(3,:));

    i = msh.LOCSTRUCT(ell,4);
    F1s_tildetheta =  (0.5*h) * sum(PLATE_WTS.*f1svals.*theta(4,:));
    F1s(i) = F1s(i) + F1s_tildetheta;

    F2s_tildetheta =  (0.5*h) * sum(PLATE_WTS.*f2svals.*theta(4,:));
    F2s(i) = F2s(i) + F2s_tildetheta;
    
    Fc(i) = Fc(i) + (0.5*h) * fcval * sum(PLATE_WTS.*theta(4,:));
end


clear F1s_tildetheta F1s_tildethetax F2s_tildetheta F2s_tildethetax x1 glbpts f1svals f2svals





%% Arrange the LHSf and LHSs as Block Matrices

LHSs = sparse(2*dofs + 1,2*dofs + 1);
LHSf = sparse(2*Mf + Mp,2*Mf + Mp);
    
%Structure LHS
%Matrics that affect w
LHSs(1:dofs, 1:dofs) = S; %11 block
LHSs(dofs+1:2*dofs, 1:dofs) = beta*Mass_s; %21 block

%Matrices that affect w_t
LHSs(1:dofs, dofs+1:2*dofs) = beta*Mass_s; %12 block
LHSs(dofs+1:2*dofs,dofs+1:2*dofs) = -Mass_s; %22 block
LHSs(2*dofs + 1, dofs + 1:2*dofs) = Es; %32 block

%Matrices that affect \tilde{c}
LHSs(1:dofs,2*dofs+1) = -Es'; %13 block


%Fluid LHS
LHSf(1:Mf,1:Mf) = 1*beta*Massf + nu*kf; %11 block
LHSf(Mf+1:2*Mf,Mf+1:2*Mf) = 1*beta*Massf + nu*kf; %22 block
LHSf(1:Mf,2*Mf+1:2*Mf+Mp) = Bx; %13 block
LHSf(Mf+1:2*Mf,2*Mf+1:2*Mf+Mp) = By; %23 block
LHSf(2*Mf+1:2*Mf+Mp,1:Mf) = -1*transpose(Bx); %31 block
LHSf(2*Mf+1:2*Mf+Mp,Mf+1:2*Mf) = -1*transpose(By); %32 block



%% Loops through plate nodes to ensure clamped boundary conditions are satisfied

%Plate displacement
LHSs(J0s,:) = 0;
LHSs(1:2*dofs,J0s) = 0;
LHSs(J0s + Ms,:) = 0;
LHSs(1:2*dofs,J0s + Ms) = 0;

Mass_s(J0s,:) = 0;
Mass_s(:,J0s) = 0;
Mass_s(J0s + Ms,:) = 0;
Mass_s(:,J0s + Ms) = 0;

%Set plate forces to zero at corners of Omega
F1s(J0s) = 0;
F1s(J0s + Ms) = 0;

F2s(J0s) = 0;
F2s(J0s + Ms) = 0;

Fc(J0s) = 0;
Fc(J0s + Ms) = 0;


%Plate velocity
LHSs(J0s + dofs,:) = 0;
LHSs(:,J0s + dofs) = 0;
LHSs(J0s + Ms + dofs,:) = 0;
LHSs(:,J0s + Ms + dofs) = 0;


    LHSs(J0s,J0s) = eye(M0s);
    LHSs(J0s + Ms, J0s + Ms) = eye(M0s);
    LHSs(J0s + dofs, J0s + dofs) = eye(M0s);
    LHSs(J0s + dofs + Ms, J0s + dofs + Ms) = eye(M0s);

%% Applies penalty to fluid nodes to ensure Dirichlet boundary conditions are satisfied

   LHSf(J0f,J0f) = penalty * eye(M0f);
   LHSf(J0f + Mf,J0f + Mf) = penalty * eye(M0f);
    
   
%% Solving for plate
RHSs = zeros(2*dofs+1,1);

dLs = decomposition(LHSs);


    RHSs(1:dofs) = F1s + Fc;
    RHSs(dofs+1:2*dofs) = F2s;
    RHSs(2*dofs+1) = 0;
    

    coefficients_s = full(dLs\RHSs);
    wn = coefficients_s(1:dofs);
    wntilde = coefficients_s(dofs+1:2*dofs);
    cn = coefficients_s(2*dofs+1);
    
gamma2 = zeros(M0f,1);
gamma2(Omega_tilde) = wntilde(msh.HERMITE);
gamma2(msh.STRUCT3(:,3)) = wntilde(msh.LOCSTRUCT(:,1:4))*thetamid(1:4) + 0.5*h*wntilde(msh.LOCSTRUCT(:,1:2)+Ms)*thetamid(5:6);

%% Solving for fluid

RHSf = zeros(2*Mf + Mp,1);


dLHSf = decomposition(LHSf);
    
    

    RHSf(1:Mf) = F1f + F1p;
    RHSf(Mf+1:2*Mf) = F2f + F2p; 

   %Makes the RHS vector satisfy the Dirichlet data
   %This section handles interior nodes.
   RHSf(J0f) = 0;
   RHSf(J0f + Mf) = penalty * gamma2(J0f);
   
   %Solve for the coefficients of fluids and pressure
   coefficients_f = full(dLHSf\RHSf);
   vn1 = coefficients_f(1:Mf);
   vn2 = coefficients_f(Mf+1:2*Mf);
   pn = coefficients_f(2*Mf+1:2*Mf+Mp);
  

    
clear  coefficients_f
%% Plots the FEM and True Solutions

plot_femsoln = 0; %Plot the fem solution? 1 = yes, 0 = no.
plot_testsoln = 0; %Plot solution to test problem
plot_mesh = 0; %Plot the mesh? 1 = yes, 0 = no.
 
if plot_mesh == 1
    meshplotter2D
    meshplot1d
end

if plot_femsoln == 1
   
   if plot_testsoln == 1
    x = msh.POS(:,1);
    y = msh.POS(:,2);
    xs = 0:0.0005:1; 
    
   v1T = 6*(y.^2 - y).*(x-1).^3.*x.^3; %Final x component of velocity
   v2T = -3*(2*y.^3 - 3*y.^2).*x.^2.*(x-1).^2.*(2*x-1); %Final y component of velocity
   wtT = 3*xs.^2.*(xs-1).^2.*(2*xs-1); %Final plate velocity
   wT = -wtT; %Final plate displacement
   pT = msh.POS(Jp,1) - msh.POS(Jp,2) + 1;
   
   wtxT = 6*(2*xs - 1).*(xs - 1).^2.*xs + 6*(2*xs - 1).*(xs - 1).*xs.^2 + 6*(xs - 1).^2.*xs.^2;
   wxT = -wtxT; 
   cT = 1;
   
   %Test solution x velocity
        figure(5)
        trisurf(msh.TRIANGLES6(:,1:3), msh.POS(:,1),msh.POS(:,2), v1T);
        shading interp;
        axis('square');
        xlabel('x','fontsize',14)
        ylabel('y','fontsize',14)
 %      view(2);
         title('X Velocity');   
         
    %Difference between fem and test x velocity     
       figure(6)
        trisurf(msh.TRIANGLES6(:,1:3), msh.POS(:,1),msh.POS(:,2), abs(v1T-vn1));
        shading interp;
        axis('square');
        xlabel('x','fontsize',14)
        ylabel('y','fontsize',14)
        zlim([0, max(v1T)])
        caxis([0, max(v1T)])
%         view(2);
         title('X Velocity minus FEM X velocity');  
         
      %Test solution x velocity
        figure(8)
        trisurf(msh.TRIANGLES6(:,1:3), msh.POS(:,1),msh.POS(:,2), v2T);
        shading interp;
        axis('square');
        % axis([0,1,0,1]);
        xlabel('x','fontsize',14)
        ylabel('y','fontsize',14)
 %        view(2);
         title('Y Velocity');
  
       %Difference between fem and test y velocity  
        figure(9)
        trisurf(msh.TRIANGLES6(:,1:3), msh.POS(:,1),msh.POS(:,2), abs(v2T-vn2));
        shading interp;
        axis('square');
        % axis([0,1,0,1]);
        xlabel('x','fontsize',14)
        ylabel('y','fontsize',14)
        zlim([0 max(v2T)])
        caxis([0 max(v2T)])
  %       view(2);
         title('Y Velocity minus FEM Y velocity');
         
                 figure(11)
        trisurf(msh.TH3(:,1:3), msh.POS(Jp,1),msh.POS(Jp,2), pT);
        shading interp;
        axis('square');
        % axis([0,1,0,1]);
        xlabel('x','fontsize',14)
        ylabel('y','fontsize',14)
%       view(2);
        title('True Pressure');
        
                 
                 figure(16)
        trisurf(msh.TH3(:,1:3), msh.POS(Jp,1),msh.POS(Jp,2), abs(pT-(pn+cn)));
        shading interp;
        axis('square');
        % axis([0,1,0,1]);
        xlabel('x','fontsize',14)
        ylabel('y','fontsize',14)
        zlim([0, max(pT+cT)])
        caxis([0, max(pT+cT)])
%       view(2);
        title('True Pressure minus FEM Pressure');
         
                
    
 figure(12)
 plot(xs,wT,'LineWidth',2.0)
hold on 
scatter(msh.LOCPOS(:,1),wn(1:Ms),'filled','k')
title('Plate Displacement')
  legend('Exact solution', 'FEM solution')
hold off

figure(13)
  plot(xs,wxT,'LineWidth',2.0)
title('Plate Displacement (Hermite Points)')
hold on
scatter(msh.LOCPOS(msh.HERMITE,1),wn(Ms+1:dofs),'filled','k')
  legend('Exact solution', 'FEM solution')
hold off

figure(14)
 plot(xs,wtT,'LineWidth',2.0)
hold on 
scatter(msh.LOCPOS(:,1),wntilde(1:Ms),'filled','k')
 title('Plate Velocity')
 legend('Exact solution', 'FEM solution')
hold off

figure(15)
 plot(xs,wtxT,'LineWidth',2.0)
title('Plate Velocity (Hermite Points)')
hold on
scatter(msh.LOCPOS(msh.HERMITE,1),wntilde(Ms+1:dofs),'filled','k')
 legend('Exact solution', 'FEM solution')
hold off
  
        
   end
       %Fem x velocity
       figure(4)
        trisurf(msh.TRIANGLES6(:,1:3), msh.POS(:,1),msh.POS(:,2), vn1);
        shading interp;
        axis('square');
        xlabel('x','fontsize',14)
        ylabel('y','fontsize',14)
 %        view(2);
         title('FEM X Velocity');
        
     %fem y velocity
     figure(7)
        trisurf(msh.TRIANGLES6(:,1:3), msh.POS(:,1),msh.POS(:,2), vn2);
        shading interp;
        axis('square');
        % axis([0,1,0,1]);
        xlabel('x','fontsize',14)
        ylabel('y','fontsize',14)
  %       view(2);
         title('FEM Y Velocity');

   %pressure     
   figure(10)
        trisurf(msh.TH3(:,1:3), msh.POS(Jp,1),msh.POS(Jp,2), pn + cn);
        shading interp;
        axis('square');
        % axis([0,1,0,1]);
        xlabel('x','fontsize',14)
        ylabel('y','fontsize',14)
 %        view(2);
         title('FEM Pressure');


end

clear plot_femsoln plot_testsoln plot_initdata plot_mesh

%% Compute the norm
Spectral_Norm

Solnnorm(blah) = HPhinorm;
Fnorm(blah) = HPhistarnorm;

end

%% Plot the Gevrey Class

BETA = transpose(1:1:trials);

figure(16)
plot(log(BETA), log(Fnorm./Solnnorm))
title(['$\log\left(\frac{\|\Phi^*\|}{\|\Phi\|}\right) ' ...
    '\mbox{ vs }\log(|\beta|)$'], interpreter = 'latex')
xlabel('$\log(|\beta|)$', interpreter = 'latex')
ylabel('$\log\left(\frac{\|\Phi^*\|}{\|\Phi\|}\right)$', interpreter = 'latex')

GC = log(Fnorm./Solnnorm)./log(BETA);
   