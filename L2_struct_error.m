% This program computes the L2-norm error of the plate displacement



%% Error computation

L2plterror = 0;
 
 for ell = 1:1:Ls
      w_I = zeros(1, nbPlateintpts);
      x1 = msh.LOCPOS(msh.LOCSTRUCT(ell,1),1);
      glbpts = (x1+h/2)*ones(1,nbPlateintpts) + h/2*PLATE_INTPTS;
      wvals = -exp(-Tend)*3*glbpts.^2.*(glbpts-1).^2.*(2*glbpts-1);
    for r = 1:1:Ts
       i = msh.LOCSTRUCT(ell,r);
       w_I = w_I + wn(i)*theta(r, :);
    end
    for r = 1:2
       i = msh.LOCSTRUCT(ell,r);
       w_I = w_I + 0.5*h*wn(i+Ms)*theta(r+4, :);
    end
    diff = (wvals - w_I).^2;
    errorchunk = 0.5*h * sum(PLATE_WTS.*diff);
    L2plterror = L2plterror + errorchunk;
 end

 L2plterror = sqrt(L2plterror);