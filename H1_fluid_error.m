
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

phixAB = zeros(Tf,3);
phiyAB = zeros(Tf,3);

ABPTS = [0,0.5;
    0.5,0.5;
    0.5,0];

for i = 1:1:3
    phixAB(1,i) = phi1x(ABPTS(i,1),ABPTS(i,2));
    phixAB(2,i) = phi2x(ABPTS(i,1),ABPTS(i,2));
    phixAB(3,i) = phi3x(ABPTS(i,1),ABPTS(i,2));
    phixAB(4,i) = phi4x(ABPTS(i,1),ABPTS(i,2));
    phixAB(5,i) = phi5x(ABPTS(i,1),ABPTS(i,2));
    phixAB(6,i) = phi6x(ABPTS(i,1),ABPTS(i,2));
    
    phiyAB(1,i) = phi1y(ABPTS(i,1),ABPTS(i,2));
    phiyAB(2,i) = phi2y(ABPTS(i,1),ABPTS(i,2));
    phiyAB(3,i) = phi3y(ABPTS(i,1),ABPTS(i,2));
    phiyAB(4,i) = phi4y(ABPTS(i,1),ABPTS(i,2));
    phiyAB(5,i) = phi5y(ABPTS(i,1),ABPTS(i,2));
    phiyAB(6,i) = phi6y(ABPTS(i,1),ABPTS(i,2));
end

clear phi1x phi2x phi3x phi4x phi5x phi6x
clear phi1y phi2y phi3y phi4y phi5y phi6y


%% H1-error calculation

v1xtest = @(x,y) 18*exp(-Tend)*(y^2-y)*((x-1)^2*x^3 + x^2*(x-1)^3);
v1ytest = @(x,y) 6*exp(-Tend)*x^3*(x-1)^3*(2*y-1);

v2xtest = @(x,y) -6*exp(-Tend)*(2*y^3-3*y^2)*(x*(2*x-1)*(x-1)^2 + x^2*(x-1)*(2*x-1)+x^2*(x-1)^2);
v2ytest = @(x,y) -18*exp(-Tend)*(y^2-y)*((x-1)^2*x^3 + x^2*(x-1)^3); 

% Int11 = 0;
% Int21 = 0;
% Int31 = 0;
% 
% Int12 = 0;
% Int22 = 0;
% Int32 = 0;

Int1 = 0;
Int2 = 0;

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

    v1x = [v1xtest(x_4,y_4), v1xtest(x_5,y_5), v1xtest(x_6,y_6)]; 
    v1y = [v1ytest(x_4,y_4), v1ytest(x_5,y_5), v1ytest(x_6,y_6)];

    v2x = [v2xtest(x_4,y_4), v2xtest(x_5,y_5), v2xtest(x_6,y_6)]; 
    v2y = [v2ytest(x_4,y_4), v2ytest(x_5,y_5), v2ytest(x_6,y_6)];
    
    D = (x_3 - x_1)*(y_2 - y_1) - (x_2 - x_1)*(y_3 - y_1);
    
    tildex_x = (y_2 - y_1)/D;
    tildex_y = -1*(x_2 - x_1)/D;
    tildey_x = -1*(y_3 - y_1)/D;
    tildey_y = (x_3 - x_1)/D;
    
%     Int11 = Int11 + abs(D)/6*(sum(v1x.^2) + sum(v1y.^2));
%     Int12 = Int12 + abs(D)/6*(sum(v2x.^2) + sum(v2y.^2));

    v1Ix = 0;
    v1Iy = 0;

    v2Ix = 0;
    v2Iy = 0;
    
    
    for r = 1:1:Tf
        v1Ix = v1Ix + vn1(msh.TRIANGLES6(ell,r))*(tildex_x*phixAB(r,:) + tildey_x*phiyAB(r,:));
        v1Iy = v1Iy + vn1(msh.TRIANGLES6(ell,r))*(tildex_y*phixAB(r,:) + tildey_y*phiyAB(r,:));
        
        v2Ix = v2Ix + vn2(msh.TRIANGLES6(ell,r))*(tildex_x*phixAB(r,:) + tildey_x*phiyAB(r,:));
        v2Iy = v2Iy + vn2(msh.TRIANGLES6(ell,r))*(tildex_y*phixAB(r,:) + tildey_y*phiyAB(r,:));
        
    end
    
%     Int21 = Int21 - 2*abs(D)/6*(sum(v1x.*v1Ix) + sum(v1y.*v1Iy));
%     Int31 = Int31 + abs(D)/6*(sum(v1Ix.^2) + sum(v1Iy.^2));

%     Int22 = Int22 - 2*abs(D)/6*(sum(v2x.*v2Ix) + sum(v2y.*v2Iy));
%     Int32 = Int32 + abs(D)/6*(sum(v2Ix.^2) + sum(v2Iy.^2));

    Int1 = Int1 + abs(D)/6*(sum((v1x-v1Ix).^2) + sum((v1y-v1Iy).^2));
    Int2 = Int2 + abs(D)/6*(sum((v2x-v2Ix).^2) + sum((v2y-v2Iy).^2));
end

% H1ferror1 = sqrt(Int11 + Int21 + Int31);
% H1ferror2 = sqrt(Int12 + Int22 + Int32);

 H1ferror1 = sqrt(Int1);
 H1ferror2 = sqrt(Int2);

%disp(H1error);
