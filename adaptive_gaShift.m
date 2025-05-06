function sp = adaptive_gaShift(E,A,V,gaParam)

sA = V'*A*V; sE  = V'*E*V; sEig = eig(-full(sA),-full(sE));
sminr = min(real(sEig)); smaxr = max(real(sEig));smini = min(imag(sEig)); smaxi = max(imag(sEig));
smin = complex(sminr, smini); smax = complex(smaxr, smaxi);

MaxIt = gaParam.MaxIt; nPop = gaParam.npop; VarSize = gaParam.Varsize;

pC = gaParam.pC; beta = gaParam.beta; mu = gaParam.mu;

nC = round(pC*nPop/2)*2; sigma = gaParam.sigma; gamma = gaParam.gamma; dp=gaParam.dp;

individual.Pos = []; individual.fit = [];

best_sol.fit = complex(Inf,Inf);

Pop = repmat(individual, nPop, 1);

% Initiate the Search Space%
%--------------------
for k = 1:nPop
    if(length(sEig)==1)
        Pop(k).Pos = sEig;
    else
        rp = unifrnd(real(smin),real(smax),VarSize);
        ip = unifrnd(imag(smin),imag(smax),VarSize);
        Pop(k).Pos = complex(rp, ip);
    end
    Kr = max(real(Pop(k).Pos)); Ki = max(imag(Pop(k).Pos));
    K = complex(Kr, Ki); 
    uno = ones(1,size(Pop(k).Pos,2));
    
    Pop(k).fit = newpole(sort([smax;smin.';Pop(k).Pos']),Pop(k).Pos',K*uno',dp);
    
    if Pop(k).fit<best_sol.fit
        best_sol = Pop(k);
    end
  
end


bestfit = nan(MaxIt,1);

% Main Loop%
%--------------------

for j = 1:MaxIt
    
    c = [Pop.fit];
    avgc = mean(c);
    
    if avgc~=0
        c = c/avgc;
    end
    
    probs = exp(-beta*c);
    
    Popc = repmat(individual,nC/2, 2);
    
    % Crossover Loop %
    %--------------------
    for t=1:nC/2
        
        [p1, p2] = RWS(Pop, probs); 
        
        [Popc(t,1).Pos, Popc(t,2).Pos] = uniCross(Pop(p1).Pos,Pop(p2).Pos,gamma);
        
    end
    
    Popc = Popc(:);
    
    % Mutation Loop %
    %--------------------
    
    for l = 1:nC
        
        Popc(l).Pos = Mutation(Popc(l).Pos, mu, sigma);
        T = max(best_sol.Pos);
        dno = ones(1,size(Popc(l).Pos,1));
        
        Popc(l).Pos = max(Popc(l).Pos, smin); Popc(l).Pos = min(Popc(l).Pos, smax);
        
        Popc(l).fit = newpole(sort([smax;smin.';Popc(l).Pos']),Popc(l).Pos',T*dno',dp);
        
        if Popc(l).fit<best_sol.fit
            best_sol = Popc(l);
        end
        
    end
    
    Pop = [Pop;Popc]; [~, so]= sort([Pop.fit]);
    Pop = Pop(so);
    
    Pop = Pop(1:nPop); 
    
    bestfit(j) = best_sol.fit;
    
    fprintf('GA- Iteration No:%d, Best Fit Shift=%d\n',j,bestfit(j));
        
end

sp = [Pop.fit];

return
