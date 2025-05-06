function [Sol,nrmrestol,tol_time]= rksm_ga(E,A,B,iter,r,tol,tolY,gaParam)

tic;
m = iter; p = size(B,2); I = speye(p); O=0*I;
[V,irr]=qr(B,0);
rr=inv(irr);nrmb = norm(inv(rr),'fro')^2;

beta = V'*B;

fprintf(' no its backward error\n')

VV= zeros(r,p*(m+2)); VV(1:r,1:p)=V; H = zeros(p*(m+2),p*(m+1));
nrmrestol=[]; nrma = norm(A,'fro');

if (norm(E-speye(r),'fro')>tolY)
    condestE = condest(E);
    singE = condestE/norm(E,'fro');
else
    singE = 1;
end

Atil = V'*A*V; Etil = V'*A*V;

i = 0; h=0;

Gsp = adaptive_gaShift(E,A,V,gaParam); Gsp = Gsp(:); Len = length(Gsp);

while i<m
     i = i+1;
    if h<Len
        h=h+1;
    else
        Gsp1 = adaptive_gaShift(E,A,V,gaParam);
        Gsp = [Gsp;Gsp1(:)]; 
        Len = length(Gsp);
        h= h+1;
    end
wrk = (A- Gsp(h)*E)\V;

    % -------------------------
    jms = (i- 1)*p+ 1; j1s = (i+1)*p; js=i*p; js1= js+1;
    
    for it = 1:2
        for kk = 1:i
            k1 = (kk-1)*p+1;
            k2 = kk*p;
            gamma = VV(1:r,k1:k2)'*wrk; H(k1:k2,jms:js)= H(k1:k2,jms:js)+ gamma;
            wrk = wrk- VV(:,k1:k2)*gamma;
        end
    end
    
    [V, H(js1:j1s,jms:js)] = qr(wrk,0);
    ih1 = i+ 1; ih = i;
    newAv = A*V;
    g = (VV(1:r,1:js)'*A*V);
    
    Ttil = Etil\Atil; btil = Etil\beta;
    
    Y = lyap(Ttil,btil*btil'); nrmx = norm(Y,'fro');
    
    % Compute Residual
    % -------------------------
    
    u1 = newAv- VV(1:r,1:js)*g;
    d=-VV(1:r,1:js)*(Y*(H(1:ih*p,1:ih*p)'\[sparse(p*(ih-1),p);I])*H(p*ih+1:p*ih1,p*ih-p+1:p*ih)');
    U=[-V*Gsp(h),  d u1 ];
    rr=qr(full(U),0); rr = triu(rr(1:size(rr,2),:));
    nrmres=norm(rr*sparse([O I O; I O I; O I O ])*rr','fro')/(nrmb+singE*nrma*nrmx);
    nrmrestol=[nrmrestol,nrmres];
    disp([i,nrmres])
    if (nrmres<tol)
        break
    end
    
    g1 =g; g2 = V'*A*VV(1:r,1:js); g3 = V'*A*V;
    ge1 = VV(1:r,1:js)'*E*V; ge2 = V'*E*VV(1:r,1:js); ge3 = V'*E*V;
    
    Atil = [Atil g1;g2 g3]; Etil = [Etil ge1;ge2 ge3];
    beta1 = V'*B; beta= [beta;beta1];
    VV(1:r,js+1:j1s)= V;
    
end

% Reduced Rank of the Solution
%---------------------------------------
[uY, sY] = eig(Y);
[sY,id]= sort(diag(sY)); sY = flipud(sY);
uY = uY(:,id(end:-1:1));
is = sum(abs(sY));
Y0 = uY*diag(sqrt(sY));
Sol = VV(:,1:size(Y0,1))*Y0;
sz = size(Sol);
tol_time = toc;

fprintf('Space dim %d Solution rank %d time %d\n',j1s,is,tol_time)
fprintf('GA Shift(RKSM): Size of the low rank factor %d\n',sz);

return


