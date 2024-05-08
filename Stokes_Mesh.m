%This program inputs a square mesh and outputs the Taylor Hood mesh
 
% clear
% close all

%% Load in the square mesh
%  Square41
%    Square145
%    Square545
     Square2113

%These meshes throw memory errors
%   Square8321
%    Square33025

%% Define the relevant mesh parameters

Mf = msh.nbNod;
Lf = length(msh.TRIANGLES6(:,1));
Tf = length(msh.TRIANGLES6(1,:))-1;

J0f = msh.LINES3(:,1:3);
J0f = unique(J0f);
M0f = length(J0f);

Jf = transpose(setdiff(1:Mf,J0f));
Nf = length(Jf);

%% Set up the pressure mesh.

Jp = unique(msh.TRIANGLES6(:,1:3)); %Index set to translate between pressure nodes and fluid nodes
Mp = length(Jp);

msh.TH3 = zeros(Lf,4);

for ell = 1:1:Lf
    for r = 1:1:3
        msh.TH3(ell,r) = find(Jp == msh.TRIANGLES6(ell,r));
    end
end

