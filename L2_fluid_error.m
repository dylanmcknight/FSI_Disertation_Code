
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



%% H1-error calculation

v1test = @(x,y) 6*exp(-Tend)*(y.^2 - y).*(x-1).^3.*x.^3;

v2test = @(x,y) -3*exp(-Tend)*(2*y.^3 - 3*y.^2).*x.^2.*(x-1).^2.*(2*x-1);

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

    v1 = [v1test(x_4,y_4), v1test(x_5,y_5), v1test(x_6,y_6)]; 
    v2 = [v2test(x_4,y_4), v2test(x_5,y_5), v2test(x_6,y_6)]; 
    
    D = (x_3 - x_1)*(y_2 - y_1) - (x_2 - x_1)*(y_3 - y_1);
    
%     Int1 = Int1 + abs(D)/6*sum(v1.^2);
%     Int2 = Int2 + abs(D)/6*sum(v2.^2);

    v1I = 0;
    v2I = 0;
    
    
    for r = 1:1:Tf
        v1I = v1I + vn1(msh.TRIANGLES6(ell,r))*phiAB(r,:);    
        v2I = v2I + vn2(msh.TRIANGLES6(ell,r))*phiAB(r,:); 
    end
    
    Int1 = Int1 + abs(D)/6*sum((v1 - v1I).^2);
    Int2 = Int2 + abs(D)/6*sum((v2 - v2I).^2);

end

L2ferror1 = sqrt(Int1);
L2ferror2 = sqrt(Int2);

%disp(H1error);

clear Int1 Int2
