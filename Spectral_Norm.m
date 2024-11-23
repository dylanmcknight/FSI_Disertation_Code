%This Program computes the error of the spectral problem i\beta I-\scrA
%\Phi = \Phi*


phi1 = @(x,y) (1-x-y)*(1-2*x-2*y);
phi2 = @(x,y) y*(2*y-1);
phi3 = @(x,y) x*(2*x-1);
phi4 = @(x,y) 4*(1-x-y)*y;
phi5 = @(x,y) 4*x*y;
phi6 = @(x,y) 4*(1-x-y)*x;

phiAB = zeros(Tf,3);

ABPTS = [0,0.5;
    0.5,0.5;
    0.5,0];

for i = 1:1:3
    phiAB(1,i) = phi1(ABPTS(i,1),ABPTS(i,2));
    phiAB(2,i) = phi2(ABPTS(i,1),ABPTS(i,2));
    phiAB(3,i) = phi3(ABPTS(i,1),ABPTS(i,2));
    phiAB(4,i) = phi4(ABPTS(i,1),ABPTS(i,2));
    phiAB(5,i) = phi5(ABPTS(i,1),ABPTS(i,2));
    phiAB(6,i) = phi6(ABPTS(i,1),ABPTS(i,2));
end

clear phi1 phi2 phi3 phi4 phi5 phi6



%% L2 Fluid Norm calculation
   

Int1 = 0;
Int2 = 0;

F1fnorm = 0;
F2fnorm = 0;

for ell = 1:Lf
    x_1 = msh.POS(msh.TRIANGLES6(ell,1),1);
    y_1 = msh.POS(msh.TRIANGLES6(ell,1),2);
    
    x_2 = msh.POS(msh.TRIANGLES6(ell,2),1);
    y_2 = msh.POS(msh.TRIANGLES6(ell,2),2);
    
    x_3 = msh.POS(msh.TRIANGLES6(ell,3),1);
    y_3 = msh.POS(msh.TRIANGLES6(ell,3),2);
    
    x_4 = msh.POS(msh.TRIANGLES6(ell,4),1);
    y_4 = msh.POS(msh.TRIANGLES6(ell,4),2);
    
    x_5 = msh.POS(msh.TRIANGLES6(ell,5),1);
    y_5 = msh.POS(msh.TRIANGLES6(ell,5),2);
    
    x_6 = msh.POS(msh.TRIANGLES6(ell,6),1);
    y_6 = msh.POS(msh.TRIANGLES6(ell,6),2);

    D = (x_3 - x_1)*(y_2 - y_1) - (x_2 - x_1)*(y_3 - y_1);
    
%     Int1 = Int1 + abs(D)/6*sum(v1.^2);
%     Int2 = Int2 + abs(D)/6*sum(v2.^2);

  glbpts = [x_3-x_1, x_2-x_1; y_3-y_1, y_2-y_1]*transpose(FLUID_INTPTS) + [x_1;y_1];
    glbpts = transpose(glbpts);
      

    %v2 Nonzero Dirichlet B.C.
     x = glbpts(:,1);
     y = glbpts(:,2);
     f1fvals = 6*beta*1*(y.^2 - y).*(x-1).^3.*x.^3 + nu*(-12*(x - 1).^3.*x.^3 - 36*(y.^2 - y).*(x - 1).^3.*x - 108*(y.^2 - y).*(x - 1).^2.*x.^2 - 36*(y.^2 - y).*(x - 1).*x.^3) + 1;
     f2fvals = -3*beta*1*(2*y.^3 - 3*y.^2).*x.^2.*(x-1).^2.*(2*x-1) + nu*(18*(x - 1).^3.*x.^2.*(2*y - 1) + 18*(x - 1).^2.*x.^3.*(2*y - 1) + 6*(2*y.^3 - 3*y.^2).*(x - 1).^3 + 54*(2*y.^3 - 3*y.^2).*(x - 1).^2.*x + 54*(2*y.^3 - 3*y.^2).*(x - 1).*x.^2 + 6*(2*y.^3 - 3*y.^2).*x.^3) - 1;

    v1I = 0;
    v2I = 0;
        
    for r = 1:1:Tf
        i = msh.TRIANGLES6(ell,r);
        v1I = v1I + vn1(i)*phiAB(r,:);    
        v2I = v2I + vn2(i)*phiAB(r,:);

        F1fnorm = F1fnorm + abs(D)/2 * sum(FLUID_WTS.*phi(r,:)*(abs(f1fvals)).^2);
        F2fnorm = F2fnorm + abs(D)/2 * sum(FLUID_WTS.*phi(r,:)*(abs(f2fvals)).^2);
    end
    
    Int1 = Int1 + abs(D)/6*sum(abs(v1I).^2);
    Int2 = Int2 + abs(D)/6*sum(abs(v2I).^2);

end

L2fnorm1 = sqrt(Int1);
L2fnorm2 = sqrt(Int2);

F1fnorm = sqrt(F1fnorm);
F2fnorm = sqrt(F2fnorm);

clear Int1 Int2

%% L2 Plate Velocity Error computation

L2pltvelnorm = 0;

Fpltvelnorm = 0;
 
 for ell = 1:1:Ls
      wtilde_I = zeros(1, nbPlateintpts);

      x1 = msh.LOCPOS(msh.LOCSTRUCT(ell,1),1);
      glbpts = (x1+h/2)*ones(1,nbPlateintpts) + h/2*PLATE_INTPTS; %transformation from [-1,1] to [x1,x2] is x = h/2(tildex) + (x1+h/2)
        
      %Forcing of plate velocity equation
      f1svals = 3*beta*1*glbpts.^2.*(glbpts-1).^2.*(2*glbpts-1) - 720*glbpts + 360 - 1;

    for r = 1:1:Ts
       i = msh.LOCSTRUCT(ell,r);
       wtilde_I = wtilde_I + wn(i)*theta(r, :);
    end
    for r = 1:2
       i = msh.LOCSTRUCT(ell,r);
       wtilde_I = wtilde_I + 0.5*h*wntilde(i+Ms)*theta(r+4, :);

       F1s_tildetheta =  (0.5*h) * sum(PLATE_WTS.*abs(f1svals).^2.*theta(r,:));
       F1s_tildethetax =  (0.5*h)^2 * sum(PLATE_WTS.*abs(f1svals).^2.*theta(r+4,:));

       Fpltvelnorm = Fpltvelnorm  + F1s_tildetheta;
       Fpltvelnorm = Fpltvelnorm + F1s_tildethetax;
    end
    diff = (abs(wtilde_I)).^2;
    errorchunk = 0.5*h * sum(PLATE_WTS.*diff);
    L2pltvelnorm = L2pltvelnorm + errorchunk;

    F1s_tildetheta =  (0.5*h) * sum(PLATE_WTS.*abs(f1svals).^2.*theta(3,:));
    Fpltvelnorm = Fpltvelnorm + F1s_tildetheta;

    F1s_tildetheta =  (0.5*h) * sum(PLATE_WTS.*abs(f1svals).^2.*theta(4,:));
    Fpltvelnorm = Fpltvelnorm + F1s_tildetheta;
 end

 L2pltvelnorm = sqrt(L2pltvelnorm);
 Fpltvelnorm = sqrt(Fpltvelnorm);

 %% H2 Plate displacement norm computation

H2pltnorm = 0;
Fpltnorm = 0;
 
 for ell = 1:1:Ls
      w_xxI = zeros(1, nbPlateintpts);

      x1 = msh.LOCPOS(msh.LOCSTRUCT(ell,1),1);
      glbpts = (x1+h/2)*ones(1,nbPlateintpts) + h/2*PLATE_INTPTS; %transformation from [-1,1] to [x1,x2] is x = h/2(tildex) + (x1+h/2)
        
      %2nd derivative of forcing of plate displacement equation
      f2sxxvals = -3*(beta + 1)*(40*glbpts.^3 - 60*glbpts.^2 + 24*glbpts -2);


    for r = 1:1:Ts
       i = msh.LOCSTRUCT(ell,r);
       w_xxI = w_xxI + 1/(0.5*h)^2*wn(i)*thetaxx(r, :);
    end
    for r = 1:2
       i = msh.LOCSTRUCT(ell,r);
       w_xxI = w_xxI + 1/(0.5*h)*wn(i+Ms)*thetaxx(r+4, :);

       F2s_tildetheta =  (0.5*h) * sum(PLATE_WTS.*abs(f2sxxvals).^2.*theta(r,:));
       F2s_tildethetax =  (0.5*h)^2 * sum(PLATE_WTS.*abs(f2sxxvals).^2.*theta(r+4,:));

       Fpltnorm = Fpltnorm  + F2s_tildetheta;
       Fpltnorm = Fpltnorm + F2s_tildethetax;
    end
    diff = (abs(w_xxI)).^2;
    errorchunk = 0.5*h * sum(PLATE_WTS.*diff);
    H2pltnorm = H2pltnorm + errorchunk;

    F2s_tildetheta =  (0.5*h) * sum(PLATE_WTS.*abs(f2sxxvals).^2.*theta(3,:));
    Fpltvelnorm = Fpltvelnorm + F2s_tildetheta;

    F2s_tildetheta =  (0.5*h) * sum(PLATE_WTS.*abs(f2sxxvals).^2.*theta(4,:));
    Fpltvelnorm = Fpltvelnorm + F2s_tildetheta;
 end

 H2pltnorm = sqrt(H2pltnorm);

 Fpltnorm = sqrt(Fpltnorm);

 %% Compute the total error

HPhinorm = sqrt(H2pltnorm^2 + L2pltvelnorm^2 + L2fnorm1^2 + L2fnorm2^2);
HPhistarnorm = sqrt(Fpltvelnorm^2 + F1fnorm^2 + F2fnorm^2 + Fpltnorm^2);