% This program computes the H2-seminorm error of the plate displacement



%% Error computation

H2plterror = 0;
 
 for ell = 1:1:Ls
      w_xxI = zeros(1, nbPlateintpts);
      x1 = msh.LOCPOS(msh.LOCSTRUCT(ell,1),1);
      glbpts = (x1+h/2)*ones(1,nbPlateintpts) + h/2*PLATE_INTPTS;
      wxxvals = -6*exp(-Tend)*(20*glbpts.^3-30*glbpts.^2+12*glbpts-ones(1,nbPlateintpts));
    for r = 1:1:Ts
       i = msh.LOCSTRUCT(ell,r);
       w_xxI = w_xxI + 1/(0.5*h)^2*wn(i)*thetaxx(r, :);
    end
    for r = 1:2
       i = msh.LOCSTRUCT(ell,r);
       w_xxI = w_xxI + 1/(0.5*h)*wn(i+Ms)*thetaxx(r+4, :);
    end
    diff = (wxxvals - w_xxI).^2;
    errorchunk = 0.5*h * sum(PLATE_WTS.*diff);
    H2plterror = H2plterror + errorchunk;
 end

 H2plterror = sqrt(H2plterror);