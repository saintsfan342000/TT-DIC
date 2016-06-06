%F=[    0.9972   -0.1926 ;   -0.0194    1.2652];

F = [0.689334 -2.362283 ; -0.015362 1.088183];

%From Dr K
phi = (1/2)*(transpose(F)*F-eye(2));
prinvals = eig(phi);
d1 = max(prinvals);
d2 = min(prinvals);
L1 = sqrt( 2*d1 + 1 );
L2 = sqrt( 2*d2 + 1 );
[L2;L1]
%1/(L1*L2)

%FROM CODE, mod
U = sqrtm(transpose(F)*F); 
diagU= eig(U)    % <------
%1/(diagU(1)*diagU(2))

%From Dr K, mod
phi = (1/2)*(sqrtm(transpose(F)*F)-eye(2,2));
prinvals = eig(phi);
d1 = max(prinvals);
d2 = min(prinvals);
L1 = sqrt(2*d1+1);
L2 = sqrt(2*d2+1);
[L2;L1]

%FROM CODE
U = transpose(F)*F; %stretching tensor???
diagU= eig(U); %principal stretch
diagU = sqrt(diagU)     % <------


eig(F);
sqrt(eig(F))
eig(F'*F)
sqrt(eig(F'*F))  % <------