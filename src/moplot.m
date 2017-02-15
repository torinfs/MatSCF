function moplot( atoms, xyz, out, iMO, level )
%Plots a isosurface of the spatial molecular orbital specified by iMO
%
%
%
%   Inputs:
%       atoms --    1xK vector of element numbers; e.g. [6 8] for CO
%       xyz --      Kx3 coordinate matrix of nuclear positions in Angstrom
%       out --      output structure of mocalc()
%       iMO --      index specifying which MO to plot, ranked in ascending
%                     energy order.
%       level --    Contour level for the isosurface, in Angstrom^(-3/2)
%

%Find largest distance between points
maxdiff = 0;
for i = 1:numel(atoms)
    for j = i:numel(atoms)
        xdiff = abs(xyz(i,1) - xyz(j,1));
        ydiff = abs(xyz(i,2) - xyz(j,2));
        zdiff = abs(xyz(i,3) - xyz(j,3));
        maxdiff = max([xdiff, ydiff, zdiff, maxdiff]);
    end
end

%This is still not quite a satisfactory way to make my box.
l = linspace(-maxdiff*10,maxdiff*10,251);
[x, y, z] = meshgrid(l);
a0 = 0.52917721067; %A/a0

aulevel = level*a0^(3/2); %A^(-3/2)*(A/a0)^(3/2) = a0^(-3/2)

MO = zeros(size(x));

for i = 1:numel(out.basis)
    
    primsum = zeros(size(MO));
    xminA = x - out.basis(i).A(1);
    yminA = y - out.basis(i).A(2);
    zminA = z - out.basis(i).A(3);

    for k = 1:numel(out.basis(i).d)
        primsum = primsum + out.basis(i).d(k)*out.basis(i).N(k)*...
            exp(-out.basis(i).alpha(k)*(xminA.^2 + yminA.^2 +zminA.^2));
    end
    
    base = xminA.^(out.basis(i).a(1)).*yminA.^(out.basis(i).a(2)).*...
        zminA.^(out.basis(i).a(3)).*primsum;
    
    MO = MO + out.C(i, iMO)*base;
end

%Convert x, y, z to angstrom
x = x*a0; y = y*a0; z = z*a0;
figure;

surfpos = patch(isosurface(x, y, z, MO,  aulevel));
surfpos.FaceColor = 'red';
surfpos.EdgeColor = 'none';

hold on;

surfneg = patch(isosurface(x, y, z, MO, -aulevel));
surfneg.FaceColor = 'blue';
surfneg.EdgeColor = 'none';

[xatom, yatom, zatom] = sphere;
for j = 1:numel(atoms)
    atomsurf = surf(0.2*xatom+(xyz(j,1)/a0),0.2*yatom+(xyz(j,2)/a0),...
        0.2*zatom+(xyz(j,3)/a0));
    atomsurf.EdgeColor = 'none';
    atomsurf.FaceColor = [0.5, 0.5, 0.5];
end

xlabel('x (Angstrom)'); ylabel('y (Angstrom)'); zlabel('z (Angstrom)');
view(3); camlight left; axis equal;


end

