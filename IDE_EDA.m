function [fitness, p, A, fA, nA, Asize, MF, MCR, iM, S_CR, S_F, Chy, iChy, S_df, FES] = ...
    IDE_EDA(fhd, p, fitness, A, fA, nA, Asize, MF, MCR, iM, ...
                NP, D, S_CR, S_F, Chy, iChy, S_df, H, lb1, ub1, fnum, FES, FESmax, doEDA)

% Bounds
if isscalar(lb1), lb = lb1.*ones(D,1); else, lb = lb1(:); end
if isscalar(ub1), ub = ub1.*ones(D,1); else, ub = ub1(:); end

% Sort population by fitness (optional but consistent with original IDE-EDA)
[fitness, idx] = sort(fitness(:));
p = p(idx,:);                       % NP x D
X = p.';                            % D x NP

% Allocate
V = X;
U = X;

% --- RSP parameters ---
k = 3;
Ri = 1:NP;
Rank = k*(NP - Ri) + 1;
Pr = Rank ./ sum(Rank);

p_rate = 0.085 + 0.085 * (FES / FESmax);
pnum   = max(2, round(p_rate * NP));
pbest  = randi(pnum, 1, NP);

% Historical memory terminal
MCR(H) = 0.9;  MF(H) = 0.9;

rMem = floor(1 + H * rand(1, NP));
CR = MCR(rMem)' + 0.1*randn(1, NP);
CR((CR < 0) | (MCR(rMem)' == -1)) = 0;
CR(CR > 1) = 1;

F = zeros(1, NP);
for i = 1:NP
    while F(i) <= 0
        F(i) = MF(rMem(i)) + Chy(iChy);
        iChy = mod(iChy, numel(Chy)) + 1;
    end
end
F(F > 1) = 1;

% jSO tweaks
if FES < 0.6*FESmax, F(F>0.7) = 0.7; end
if FES < 0.2*FESmax
    FW = 0.7*F;
elseif FES < 0.4*FESmax
    FW = 0.8*F;
else
    FW = 1.2*F;
end
if FES < 0.25*FESmax
    CR(CR<0.7) = 0.7;
elseif FES < 0.5*FESmax
    CR(CR<0.6) = 0.6;
end

% Mutation indices
r1 = randsample(NP, NP, true, Pr);
r2 = randsample(NP, NP, true, Pr);

% --- DE mutation + crossover (for all i) ---
for i = 1:NP
    % ensure r1 != i
    tries=0;
    while r1(i)==i && tries<10
        r1(i)=randsample(NP,1,true,Pr); tries=tries+1;
    end
    tries=0;
    while (r2(i)==i || r2(i)==r1(i)) && tries<10
        r2(i)=randsample(NP,1,true,Pr); tries=tries+1;
    end

    useA = (nA > 0) && (rand < (nA/(nA+NP)));

    if useA
        aidx = randi(nA);
        V(:,i) = X(:,i) + FW(i).*(X(:,pbest(i)) - X(:,i)) + F(i).*(X(:,r1(i)) - A(:,aidx));
    else
        V(:,i) = X(:,i) + FW(i).*(X(:,pbest(i)) - X(:,i)) + F(i).*(X(:,r1(i)) - X(:,r2(i)));
    end

    % bounds (midpoint repair)
    V(:,i) = max(min(V(:,i), ub), lb);
    low = V(:,i) < lb; high = V(:,i) > ub;
    V(low,i)  = 0.5*(lb(low)  + X(low,i));
    V(high,i) = 0.5*(ub(high) + X(high,i));

    % binomial crossover
    jrand = randi(D);
    Ui = X(:,i);
    for j=1:D
        if rand <= CR(i) || j==jrand
            Ui(j) = V(j,i);
        end
    end
    U(:,i) = Ui;
end

% Evaluate trial population
fu = feval(fhd, U, fnum);
FES = FES + NP;

% Selection + archive/memory update
nS = 0;
for i = 1:NP
    if fu(i) < fitness(i)
        % archive push
        if nA < Asize
            A(:, nA+1) = X(:,i);
            fA(nA+1)   = fitness(i);
            nA = nA + 1;
        else
            ri = floor(1 + Asize*rand);
            A(:,ri) = X(:,i);
            fA(ri)  = fitness(i);
        end

        nS = nS + 1;
        S_CR(nS) = CR(i);
        S_F(nS)  = F(i);
        S_df(nS) = abs(fu(i) - fitness(i));

        X(:,i) = U(:,i);
        fitness(i) = fu(i);
    elseif fu(i) == fitness(i)
        X(:,i) = U(:,i);
    end
end

% Update MF/MCR
if nS > 0
    w = S_df(1:nS) ./ sum(S_df(1:nS));
    if all(S_CR(1:nS) == 0)
        MCR(iM) = -1;
    elseif MCR(iM) ~= -1
        MCR(iM) = (sum(w.*S_CR(1:nS).*S_CR(1:nS))/sum(w.*S_CR(1:nS)) + MCR(iM))/2;
    end
    MF(iM)  = (sum(w.*S_F(1:nS).*S_F(1:nS))/sum(w.*S_F(1:nS)) + MF(iM))/2;
    iM = mod(iM, H) + 1;
end

% Keep archive size sane
if ~isempty(fA)
    [fA, idxA] = sort(fA);
    A = A(:, idxA);
    Asize = NP;
    if nA > Asize
        nA = Asize;
        A = A(:,1:Asize);
        fA = fA(1:Asize);
    end
end

% --- Optional EDA part ---
if doEDA
    PApnum = max(2, ceil(0.5*NP));
    if PApnum < 2*D, PApnum = NP; end

    newN = max(1, ceil(0.9*pnum));
    sel = X(:,1:PApnum).';         % PApnum x D

    mu = mean(sel,1);
    C  = cov(sel,1);
    C  = (C + C.')/2 + 1e-12*eye(D);  % regularize

    temp_pos = mvnrnd(mu, C, newN).'; % D x newN

    % bound repair
    UB = repmat(ub,1,newN); LB = repmat(lb,1,newN);
    Ran = rand(D,newN);
    bad = (temp_pos>UB) | (temp_pos<LB);
    temp_pos(bad) = Ran(bad).*(UB(bad)-LB(bad)) + LB(bad);

    fuE = feval(fhd, temp_pos, fnum);
    FES = FES + newN;

    % merge + truncate best NP
    fuAll = [fuE(:); fitness(:)];
    XAll  = [temp_pos, X];

    [fuAll, idu] = sort(fuAll);
    X = XAll(:, idu(1:NP));
    fitness = fuAll(1:NP);
end

% Return back to row-wise
p = X.';     % NP x D
fitness = fitness(:);
end
