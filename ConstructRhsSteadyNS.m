function [r] = ConstructRhsSteadtNS(numUn, nP, nO, M, N, R, dr, dth, U, sol, Re)
    r = zeros(numUn, 1);
%
%  Boundary conditions at right: omega given
    for jrow = 2: N-1
        ijO = nO(jrow, M);
        r(ijO) = U/R + 3*U/dr;
    end  
    tic
%set Jacobian of the interior points, boundary points are zero
    for i = 2: N-1
        for j = 2:M-1
            ijO = nO(i,j);
            rij = (j-1)*dr;
            OmegaR = nO(i,j+1); OmegaL = nO(i,j-1);
            OmegaT = nO(i+1,j); OmegaB = nO(i-1,j);
            PsiR = nP(i,j+1); PsiL = nP(i,j-1);
            PsiT = nP(i+1,j); PsiB = nP(i-1,j);
            J = ((sol(OmegaR)-sol(OmegaL))*(sol(PsiT) -...
            sol(PsiB))-(sol(OmegaT)-sol(OmegaB))*(sol(PsiR)-sol(PsiL)));
            J = J/(dth*dr*4*rij);
            r(ijO) = -Re*J;
        end
    end
    toc
end
