function [Sol, nrmrestol,tol_time] = adi_ga(E,A,B,iter,tol,gaParam)

tic;

m = iter; uB = B; bnorm = norm((uB*uB'),'fro');

nrmrestol = []; i= 0; h = 0; Sol = [];

[V,~] = qr(uB,0);
Gsp = adaptive_gaShift(E,A,V,gaParam);Gsp = Gsp(:); Len = length(Gsp);

while i<m
    
    i = i+1;
    if  h<Len
        h= h+1;
    else
        [V,~] = qr(uB,0);
        Gsp1 = adaptive_gaShift(E,A,V,gaParam);
        Gsp = [Gsp;Gsp1(:)]; Len = length(Gsp);
        h = h+1;
    end
    
    X = (A+Gsp(h)*E)\uB;
    
    %%%%%%%%%%%%%%%%%%%%%%%
    if isreal(Gsp(h))
        Sol = [Sol X];
        uB = uB- (2*(Gsp(h)*E*X));
    else
        Gsp(h) = conj(Gsp(h));
        n = sqrt(2); delta = real(Gsp(h))/imag(Gsp(h));
        uB = uB- (4*real(Gsp(h))*E)*(real(X)+(delta*imag(X)));
        Sol = [Sol (n*(real(X)+(delta*imag(X)))) (n*(sqrt((delta)^2+1))*imag(X))];
        i=i+1;
    end
    %%%%%%%%%%%%%%%%%%%%%%%
    
    nrm = norm((uB*uB'),'fro')/bnorm;
    
    fprintf('Iteration No: %d, Residual: %d\n',i,nrm);
    
    nrmrestol = [nrmrestol nrm];
    
    if nrm<tol
        break
    end
    
end
tol_time = toc;
fprintf('GA Shift(ADI): Size of  the ADI Solution: %d, time: %d\n',size(Sol,2),tol_time);

return
