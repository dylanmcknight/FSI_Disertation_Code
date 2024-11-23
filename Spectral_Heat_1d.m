%Solve the spectral heat equation with 0 Dirichlet boundary data using quadratic
%basis functions.

%Is working

clear
close all

trials = 100;

L2norm = zeros(trials,1);
Fnorm = zeros(trials,1);

for m = 1:trials

%% Define parameters of the problem
% beta u - u_{xx} = f

beta = m*1i;

L = 16; %number of quadratic elements
M = 2*L+1; %number of nodes

dx = 1/L;%width of an element (characteristic length)

%% Build the mesh
msh.POS = zeros(M,1);
msh.ELEMENTS = zeros(L,3);
for ell = 1:L
    msh.ELEMENTS(ell,1) = 2*ell - 1;
    msh.ELEMENTS(ell,2) = 2*ell + 1;
    msh.ELEMENTS(ell,3) = 2*ell;

    %Position of local node 1 in element ell
    msh.POS(2*ell - 1) = dx*(ell-1);
    %Position of local node 3 in element ell
    msh.POS(2*ell) = msh.POS(2*ell - 1) + 0.5*dx;
end
msh.POS(M) = 1;

J0 = [1,M]; %Global node numbers associated to Dirichlet b.c.s 
M0 = length(J0);

J = setdiff(1:M,J0);

%% Plot the mesh
plot_mesh = 0;

if plot_mesh == 1
%This program plots 1d meshes
    %fh = figure;
    
    %Make full screen:
    %set(fh,'units','normalized','position',[0,0,1,1]); %[x y width height]
node_numbering=1;
    
    figure(1)
    plot(msh.POS(:,1),zeros(M,1));
    %axis('square')
    axis([0,1,-0.5,0.5]);
    xlabel('x')
    ylabel('y')
    title('Mesh (Local Node Numbering)')
    hold on;
    if node_numbering == 0
       for k = 1:M
            text(msh.POS(k,1),0,'');
       end
    else
      for k = 1:M
            text(msh.POS(k,1),0,num2str(k), 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Color', 'k');
      end
    end  

end

clear node_numbering plot_mesh k

%% Define the standard element (1) --- (3) --- (2)

phi1 = @(x) 0.5*(x.^2 - x);
phi2 = @(x) 0.5*(x.^2 + x);
phi3 = @(x) 1 - x.^2;

phi1x = @(x) 0.5*(2*x - 1);
phi2x = @(x) 0.5*(2*x + 1);
phi3x = @(x) -2*x;

INTPTS = [-sqrt(3/5), 0, sqrt(3/5)];
WTS = [5/9, 8/9, 5/9];

NBINTPTS = length(INTPTS);

phi = zeros(NBINTPTS, NBINTPTS);
phix = zeros(NBINTPTS, NBINTPTS);

for k = 1:NBINTPTS
    phi(:,k) = [phi1(INTPTS(k)); phi2(INTPTS(k)); phi3(INTPTS(k))];
    phix(:,k) = [phi1x(INTPTS(k)); phi2x(INTPTS(k)); phi3x(INTPTS(k))];
end

clear phi1 phi2 phi3 phi1x phi2x phi3x

%% Build the stiffness and mass matrices
K = zeros(M,M);
Mass = zeros(M,M);

for ell = 1:L
    for r = 1:3
        ir = msh.ELEMENTS(ell,r); 
        for s = 1:3
            ij = msh.ELEMENTS(ell,s);

            K(ir,ij) = K(ir,ij) + 2/dx*sum(WTS.*phix(r,:).*phix(s,:));
            Mass(ir,ij) = Mass(ir,ij) + dx/2*sum(WTS.*phi(r,:).*phi(s,:));
        end
    end
end

%% Compute forcing vector
%u(x) = (x-x^2)

F0 = zeros(M,1);

f = @(x) beta*(x - x.^2) + 2; %forcing

for ell = 1:L
    x1 = msh.POS(msh.ELEMENTS(ell,1));
    x2 = msh.POS(msh.ELEMENTS(ell,2));

    xpts = dx/2 *INTPTS + 0.5*(x1 + x2);

    fvals = f(xpts);

    for r = 1:3
        ir = msh.ELEMENTS(ell,r);

        F0(ir) = F0(ir) + dx/2*sum(WTS.*fvals.*phi(r,:));
    end
end


%% Change rows and columns to ensure Dirichlet B.C.s

LHS = beta*Mass + K;


LHS(J0,:) = 0;
LHS(:,J0) = 0;
LHS(J0,J0) = eye(M0,M0);

F0(J0) = 0;

%% Solving for coefficients

    RHS = F0;

    RHS(J0) = 0;

    alpha = LHS\RHS; %LHS^(-1)*RHS

%% Plot the solution
plot_soln = 0;

if plot_soln == 1

x = msh.POS;

utrue = x-x.^2;

figure(2)
plot(x,utrue)
hold on
scatter(x,alpha)
hold off
legend('true solution', 'FEM solution')

end

%% Norm computation

for ell = 1:L
    x1 = msh.POS(msh.ELEMENTS(ell,1));
    x2 = msh.POS(msh.ELEMENTS(ell,2));

    xpts = dx/2 *INTPTS + 0.5*(x1 + x2);

    fvals = f(xpts);

    u_I = 0;
    for r = 1:3
        ir = msh.ELEMENTS(ell,r);
        u_I = u_I + alpha(ir)*phi(r,:);

        Fnorm(m) = Fnorm(m) + dx/2*sum(WTS.*abs(fvals).^2.*phi(r,:));
    end
    diff = abs(u_I).^2;
    L2norm(m) = L2norm(m) + dx/2*sum(WTS.*diff);

end

    L2norm(m) = sqrt(L2norm(m));
    Fnorm(m) = sqrt(Fnorm(m));

end

BETA = transpose(1:1:trials);

figure(3)
plot(log(BETA), log(Fnorm./L2norm))
title(['$\log\left(\frac{\|\Phi^*\|}{\|\Phi\|}\right) ' ...
    '\mbox{ vs }\log(|\beta|)$'], interpreter = 'latex')
xlabel('$\log(|\beta|)$', interpreter = 'latex')
ylabel('$\log\left(\frac{\|\Phi^*\|}{\|\Phi\|}\right)$', interpreter = 'latex')

GC = log(Fnorm./L2norm)./log(BETA);