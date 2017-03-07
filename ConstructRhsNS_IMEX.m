function [r] = ConstructRhsNS_IMEX(numUn, nP, nO, M, N, R, dr, dth, dt, U, sol, solOld)
    r = zeros(numUn, 1);
%
%  Boundary conditions at right: omega given
    for jrow = 2: N-1
        ijO = nO(jrow, M);
        r(ijO) = U/R + 3*U/dr;
    end  
    
    for i = 2: N-1
        for j = 2:M-1
            ijO = nO(i,j);
            r(ijO) = 4*sol(ijO)/3;
            r(ijO) = r(ijO)-solOld(ijO)/3;
            r(ijO) = r(ijO)+2*dt*Jacob(dr, dth, i, j, nP, nO, solOld)/3;
            r(ijO) = r(ijO)-4*dt*Jacob(dr, dth, i, j, nP, nO, sol)/3;
        end
    end
     
end