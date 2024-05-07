
psi1 = @(x,y) 1-x-y;
psi2 = @(x,y) y;
psi3 = @(x,y) x;

psiAB = zeros(Tf,3);

ABPTS = [0,0.5;
    0.5,0.5;
    0.5,0];

for i = 1:1:3
    psiAB(1,i) = psi1(ABPTS(i,1),ABPTS(i,2));
    psiAB(2,i) = psi2(ABPTS(i,1),ABPTS(i,2));
    psiAB(3,i) = psi3(ABPTS(i,1),ABPTS(i,2));
end

clear psi1 psi2 psi3


%% H1-error calculation

%qtest = @(x,y) (1+exp(-Tend))*(x-y);
qtest = @(x,y) (1+exp(-Tend))*(x.^2-y.^2);
cact = 1;

Int = 0;

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

    qact = [qtest(x_4,y_4), qtest(x_5,y_5), qtest(x_6,y_6)]; 
    
    D = (x_3 - x_1)*(y_2 - y_1) - (x_2 - x_1)*(y_3 - y_1);
    

    qI = zeros(1,3);
    
    
    for r = 1:1:3
        qI = qI + pn(msh.TH3(ell,r))*(psiAB(r,:));   
    end
    
    Int = Int + abs(D)/6*sum((qact-qI).^2);

end

L2perror = sqrt(Int + abs(cact - cn)^2);

%disp(H1error);
