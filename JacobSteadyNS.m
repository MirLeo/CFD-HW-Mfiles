% Find the Jacobian
function J = Jacob(dr, dth, i, j, nP, nO, sol)

rij = (j-1)*dr;
OmegaR = nO(i,j+1); OmegaL = nO(i,j-1);
OmegaT = nO(i+1,j); OmegaB = nO(i-1,j);
PsiR = nP(i,j+1); PsiL = nP(i,j-1);
PsiT = nP(i+1,j); PsiB = nP(i-1,j);

J = ((sol(OmegaR)-sol(OmegaL))*(sol(PsiT) -...
    sol(PsiB))-(sol(OmegaT)-sol(OmegaB))*(sol(PsiR)-sol(PsiL)));
J = J/(dth*dr*4*rij);