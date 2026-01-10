function [bestever,fbest_hist]=MGSA_EM(fhd,layers,d,MAX_FES,funcid)
% ---------------------- Core settings ----------------------
Rpower = 1;                % distance exponent in force model
lb = -100; ub = 100;       % search space bounds
layersno = length(layers);
layercumsum = cumsum(layers);
sz = sum(layers);          % total population size (sum of all layers)

% Mutation/restart parameters (HGSA part)
mut_max = 0.7;
mut_min = 0.1;
Sg = sz;                   % stagnation threshold (you can lower, e.g., 10~50)

% phi controls global best attraction in winner update (dimension-based)
if d==50
    phi=0.04;
elseif d==30
    phi=0.02;
else
    phi=0.008;
end
% ---------------------- Layer containers ----------------------
% lvlp/lvlv store positions/velocities per layer (layer ¡Á index ¡Á dim)
lvlp = zeros(layersno, max(layers), d);
lvlv = zeros(layersno, max(layers), d);
fitnesslayer = zeros(layersno, max(layers), 1);
pblayer = zeros(layersno, max(layers), d);     % personal best (PB) per layer

% ---------------------- Initialize population ----------------------
p  = lb + (ub - lb) .* rand(sz, d);
v  = 0.1*(lb + (ub - lb) .* rand(sz, d));
pb = p;                                        % PB starts at current positions

% ---------------------- IDE-EDA memory (restart operator) ----------------------
Asize = sz;          % archive size
A = []; fA = [];     % archive positions / archive fitness
nA = 0;              % archive count
H = 5;               % memory length for MF/MCR (JADE-like)
MF = 0.3 * ones(H,1);
MCR = 0.8 * ones(H,1);
iM = 1;
S_CR = zeros(1,sz);
S_F  = zeros(1,sz);
S_df = zeros(1,sz);
Chy = cauchyrnd(0,0.1,sz+10);   % Cauchy noise pool for F sampling
iChy = 1;

fbest_hist = [];

% ---------------------- Initial fitness evaluation ----------------------
af = zeros(sz,1);
for i = 1:sz
    af(i) = feval(fhd, p(i,:)', funcid);
end
fitness = af;
bestever = min(fitness);
fbest_hist = [fbest_hist, bestever];

FES = sz;                             % number of function evaluations used
MAXGEN = floor(MAX_FES/sz);
gen = 1;

% Temporary buffers for rebuild (flatten layers back into population)
tp  = zeros(sz,d);  tv  = zeros(sz,d);
tpb = zeros(sz,d);  tf  = 9999999*zeros(sz,1);
tcount = zeros(1,sz);

% For merging losers/winners output ordering
layeridxs = zeros(layersno,2);
allpos = zeros(sz,d);
allv   = zeros(sz,d);

% Stagnation tracking (IMPORTANT: must be consistent with layer ordering)
counter = zeros(1,sz);
counterlayer = zeros(layersno, max(layers));

while(FES < MAX_FES)
% ============================================================
    % 1) Rank population and assign into layers (best ¡ú top layer)
    % ============================================================
    [~, pid] = sort(fitness);
    idx = pid(1:layercumsum(1));
    lvlp(1,1:length(idx),:)=p(idx,:);
    lvlv(1,1:length(idx),:)=v(idx,:);
    fitnesslayer(1,1:length(idx),:)=fitness(idx);
    pblayer(1,1:length(idx),:)=pb(idx,:);
    counterlayer(1,1:length(idx)) = counter(idx);
    for li=2:layersno
        idx = pid(layercumsum(li-1)+1:layercumsum(li));
        len = length(idx);
        lvlp(li,1:len,:)=p(idx,:);
        lvlv(li,1:len,:)=v(idx,:);
        fitnesslayer(li,1:len,:)=fitness(idx);
        pblayer(li,1:length(idx),:)=pb(idx,:);
        counterlayer(li,1:len) = counter(idx);
    end
    % ============================================================
    % 2) Mass computation from fitness (HGSA gravity-like weighting)
    % ============================================================
    fmax = max(fitness); fmin = min(fitness);
    if fmax == fmin
        M = ones(sz,1);
    else
        M = (fitness - fmax) ./ (fmin - fmax);
    end
    M = M ./ sum(M);
    % ============================================================
    % 3) Competitive learning inside each layer (pair losers/winners)
    % ============================================================
    
    ttidx=1;
    for li=layersno:-1:1
        lvlsize = layers(li);
        rlist = randperm(lvlsize);
        seprator = floor(lvlsize/2);
        % Pair indices inside this layer
        rpairs = [rlist(1:seprator); rlist((seprator+1):(2*seprator))]';
        %competitive learning
        mask = (fitnesslayer(li,rpairs(:,1),1) > fitnesslayer(li,rpairs(:,2),1))';
        losers = mask.*rpairs(:,1) + ~mask.*rpairs(:,2);
        winners = ~mask.*rpairs(:,1) + mask.*rpairs(:,2);
        % Random coefficients for velocity updates
        randco1 = rand(seprator, d);
        randco2 = rand(seprator, d);
        randco3 = rand(seprator, d);
        randco4 = rand(seprator, d);
        % Reshape (pull layer submatrices into 2D seprator¡Ád)
        lvlvlosert=reshape(lvlv(li,losers,:),[seprator,d]);
        lvlplosert=reshape(lvlp(li,losers,:),[seprator,d]);
        lvlpblosert=reshape(pblayer(li,losers,:),[seprator,d]);

        lvlvwinert=reshape(lvlv(li,winners,:),[seprator,d]);
        lvlpwinert=reshape(lvlp(li,winners,:),[seprator,d]);
        lvlpbwinert=reshape(pblayer(li,winners,:),[seprator,d]);
        % Sample top-layer individuals as "global best pool"
        toplvlsize = layers(1);
        indciestop = 1+mod(randperm(seprator),ones(1,seprator)*toplvlsize);
        gbpmat = reshape(lvlp(1,indciestop,:),[seprator d]);

        % ---------------- Force + acceleration term (HGSA) ----------------
        diff = lvlpwinert - lvlplosert;            % direction loser->winner
        R = vecnorm(diff, 2, 2);                   % distance per pair
        F = (rand(seprator,1).*M(li)) .* (diff ./ (R.^Rpower + eps));

        % Gravitational constant decay (exploration->exploitation)
        alpha = 10; G0 = 100;
        G = G0 * exp(-alpha*gen/MAXGEN);

        a = F * G;                                  % acceleration
        MR = mut_max - rand*(mut_max-mut_min)*(gen/MAXGEN);
        MR = min(mut_max, max(mut_min, MR));        % clamp mutation rate

        % ---------------- Loser update: learn from winner (+ PB or +a) -----
        if rand<MR
            lvlvlosert2 = randco1.*lvlvlosert;
            lvlvlosert2 = lvlvlosert2 + randco2.*(lvlpwinert - lvlplosert);
            lvlvlosert2=lvlvlosert2+randco4.*a;
            lvlplosert2 = lvlplosert + lvlvlosert2;
        else
            lvlvlosert2 = randco1.*lvlvlosert;
            lvlvlosert2 = lvlvlosert2 + randco2.*(lvlpwinert - lvlplosert);
            lvlvlosert2 = lvlvlosert2 + randco3.*(lvlpblosert - lvlplosert);
            lvlplosert2 = lvlplosert + lvlvlosert2;
        end
        % ---------------- Winner update: learn from upper layer ------------
        if li~=1
            upperlvlsize = layers(li-1);
            indcies = 1+mod(randperm(seprator),ones(1,seprator)*upperlvlsize);
            upperpmat = reshape(lvlp(li-1,indcies,:),[seprator d]);

            randco1 = rand(seprator, d);
            randco2 = rand(seprator, d);
            randco3 = rand(seprator, d);
            randco4 = rand(seprator, d);
            diff1 = upperpmat-lvlpwinert;          % seprator¡Ád
            R1 = vecnorm(diff1, 2, 2);                % seprator¡Á1
            F1 = (rand(seprator,1).*M(li)) .* (diff1 ./ (R1.^Rpower + eps));
            % Calculation of accelaration.
            a1=F1*G;

            if rand<MR
                % velocity and position update for winer
                lvlvwinert2 = randco1.*lvlvwinert;
                lvlvwinert2 = lvlvwinert2 + randco2.*(upperpmat- lvlpwinert);
                lvlvwinert2 = lvlvwinert2 +randco3.*(a1);
                lvlvwinert2=reshape(lvlvwinert2,[seprator,d]);
                lvlpwinert2 = lvlpwinert + lvlvwinert2;
            else
                lvlvwinert2 = randco1.*lvlvwinert;
                lvlvwinert2 = lvlvwinert2 + randco2.*(upperpmat - lvlpwinert);
                lvlvwinert2 = lvlvwinert2 + randco3.*(lvlpbwinert - lvlpwinert);
                lvlvwinert2 = lvlvwinert2 + phi*randco4.*(gbpmat- lvlpwinert);

                lvlvwinert2=reshape(lvlvwinert2,[seprator,d]);
                lvlpwinert2 = lvlpwinert + lvlvwinert2;
            end
        else
            lvlvwinert2 = lvlvwinert;
            lvlpwinert2 = lvlpwinert;
        end

        % combine positions of both loser and winner
        mergedlvlp=[lvlplosert2;lvlpwinert2];
        mergelen = size(lvlplosert2,1);

        ts1=mergelen; ts2=size(lvlpwinert2,1);
        layeridxs(li,:)=[ts1,ts2];

        % combine velocities of both loser and winner
        mergedlvlv = [lvlvlosert2; lvlvwinert2];
        allpos(ttidx:ttidx+(ts1+ts2-1),:)=mergedlvlp;
        allv(ttidx:ttidx+(ts1+ts2-1),:)=mergedlvlv;
        ttidx = ttidx+(ts1+ts2);
        llosers{li}=losers;
        wwiners{li}=winners;
    end
    % ============================================================
    % 4) Bound handling + evaluate ALL candidates (allpos)
    % ============================================================
    allpos(allpos>ub)=ub;
    allpos(allpos<lb)=lb;

    ffs = zeros(sz,1);
    for i=1:sz
        ffs(i) = feval(fhd, allpos(i,:)', funcid);
        if ffs(i)<fitness(i)
            fitness(i)=ffs(i);
            p(i,:)=allpos(i,:);
            counter(i)=0;
        else
            counter(i)=counter(i)+1;
        end
    end
    FES = FES + sz;
    if FES >= MAX_FES
        break;
    end

    % ============================================================
    % 5) Push improvements back into layer arrays (elitist update)
    % ============================================================
    ttidx = 1;
    for li = layersno:-1:1
        ts1 = layeridxs(li,1);
        ts2 = layeridxs(li,2);

        ff1 = ffs(ttidx:ttidx+ts1-1);
        posL = allpos(ttidx:ttidx+ts1-1,:);
        velL = allv(ttidx:ttidx+ts1-1,:);
        ff2 = ffs(ttidx+ts1:ttidx+ts1+ts2-1);
        posW = allpos(ttidx+ts1:ttidx+ts1+ts2-1,:);
        velW = allv(ttidx+ts1:ttidx+ts1+ts2-1,:);
        ttidx = ttidx + ts1 + ts2;

        losers  = llosers{li};
        winners = wwiners{li};

        % Losers improved?
        goodloseridx = fitnesslayer(li,losers,1)' > ff1;
        if any(goodloseridx)
            lgood = losers(goodloseridx);
            lvlp(li,lgood,:) = posL(goodloseridx,:);
            lvlv(li,lgood,:) = velL(goodloseridx,:);
            fitnesslayer(li,lgood,:) = ff1(goodloseridx);
        end

        % Winners improved?
        goodwinderidx = fitnesslayer(li,winners,1)' > ff2;
        if any(goodwinderidx)
            wgood = winners(goodwinderidx);
            lvlp(li,wgood,:) = posW(goodwinderidx,:);
            lvlv(li,wgood,:) = velW(goodwinderidx,:);
            fitnesslayer(li,wgood,:) = ff2(goodwinderidx);
        end

        % Stagnation counters (layer-consistent)
        counterlayer(li, losers)  = counterlayer(li, losers)  + 1;
        counterlayer(li, winners) = counterlayer(li, winners) + 1;
        counterlayer(li, losers(goodloseridx))   = 0;
        counterlayer(li, winners(goodwinderidx)) = 0;
    end

    % ============================================================
    % 6) Rebuild population arrays from layers
    % ============================================================
    for li=1:layersno
        if li==1
            leftb=1;
            rightb=layercumsum(li);
        else
            leftb=layercumsum(li-1)+1;
            rightb = layercumsum(li);
        end

        tp(leftb:rightb,:) = reshape(lvlp(li,1:layers(li),:),[layers(li) d]);
        tv(leftb:rightb,:) = reshape(lvlv(li,1:layers(li),:),[layers(li) d]);
        tpb(leftb:rightb,:) = reshape(pblayer(li,1:layers(li),:),[layers(li) d]);
        tf(leftb:rightb,:)=reshape(fitnesslayer(li,1:layers(li)),[layers(li) 1]);
        tcount(leftb:rightb) = counterlayer(li,1:layers(li));
    end

   % ============================================================
    % 7) Restart trigger: call IDE-EDA only when stagnation is high
    % ============================================================
    p = tp;  v = tv;  pb = tpb;  fitness = tf;
    counter=tcount;
    stagRate = mean(counter > Sg);
    if stagRate > 0.5   % e.g., >50% stagnated
        doEDA = true;  % you can set false to run only DE part (cheaper)

        [fitness, p, A, fA, nA, Asize, MF, MCR, iM, S_CR, S_F, Chy, iChy, S_df, FES] = ...
            IDE_EDA(fhd, p, fitness, A, fA, nA, Asize, MF, MCR, iM, ...
            sz, d, S_CR, S_F, Chy, iChy, S_df, H, lb, ub, funcid, FES, MAX_FES, doEDA);
        pb = p;
        v  = zeros(sz,d);
        counter(:) = 0;
    end


    % ============================================================
    % 8) Logging best solution
    % ============================================================
    bestever = min(bestever, min(fitness));
    fbest_hist = [fbest_hist, bestever];
    gen = gen + 1;
end
end

function r= cauchyrnd(varargin)
% % % % %Generate random numbers from the Cauchy distribution, r= a + b*tan(pi*(rand(n)-0.5)).
% % % % Chy = cauchyrnd(0, 0.1, SearchAgents_no + 10);%(SearchAgents_no + 10)*(SearchAgents_no + 10) rand numbers
% USAGE:       r= cauchyrnd(a, b, n, ...)
% Generate random numbers from the Cauchy distribution, r= a + b*tan(pi*(rand(n)-0.5)).
% ARGUMENTS:
% a (default value: 0.0) must be scalars or size(x).
% b (b>0, default value: 1.0) must be scalars or size(x).
% n and onwards (default value: 1) specifies the dimension of the output.
% EXAMPLE:
% r= cauchyrnd(0, 1, 10); % A 10 by 10 array of random values, Cauchy distributed.
% SEE ALSO:    cauchycdf, cauchyfit, cauchyinv, cauchypdf.
% Copyright (C) Peder Axensten <peder at axensten dot se>
% HISTORY:
% Version 1.0, 2006-07-10.
% Version 1.1, 2006-07-26.
% - Added cauchyfit to the cauchy package.
% Version 1.2, 2006-07-31:
% - cauchyinv(0, ...) returned a large negative number but should be -Inf.
% - Size comparison in argument check didn't work.
% - Various other improvements to check list.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Default values
a=	0.0;
b=	1.0;
n=	1;
% Check the arguments
if(nargin >= 1)
    a=	varargin{1};
    if(nargin >= 2)
        b=			varargin{2};
        b(b <= 0)=	NaN;	% Make NaN of out of range values.
        if(nargin >= 3)
            n=	[varargin{3:end}];
        end
    end
end
% Generate
r=	cauchyinv(rand(n), a, b);
end

function x= cauchyinv(p, varargin)
% USAGE:       x= cauchyinv(p, a, b)
% Inverse of the Cauchy cumulative distribution function (cdf), x= a + b*tan(pi*(p-0.5)).
% ARGUMENTS:
% p (0<=p<=1) might be of any dimension.
% a (default value: 0.0) must be scalars or size(p).
% b (b>0, default value: 1.0) must be scalars or size(p).
% EXAMPLE:
% p= 0:0.01:1;
% plot(cauchyinv(p), p);
% SEE ALSO:    cauchycdf, cauchyfit, cauchypdf, cauchyrnd.
% Copyright (C) Peder Axensten <peder at axensten dot se>
% HISTORY:
% Version 1.0, 2006-07-10.
% Version 1.1, 2006-07-26.
% - Added cauchyfit to the cauchy package.
% Version 1.2, 2006-07-31:
% - cauchyinv(0, ...) returned a large negative number but should be -Inf.
% - Size comparison in argument check didn't work.
% - Various other improvements to check list.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Default values
a=	0.0;
b=	1.0;
% Check the arguments
if(nargin >= 2)
    a=	varargin{1};
    if(nargin == 3)
        b=			varargin{2};
        b(b <= 0)=	NaN;	% Make NaN of out of range values.
    end
end
if((nargin < 1) || (nargin > 3))
    error('At least one argument, at most three!');
end
p(p < 0 | 1 < p)=	NaN;
% Calculate
x=			a + b.*tan(pi*(p-0.5));
% Extreme values.
if(numel(p) == 1), 	p= repmat(p, size(x));		end
x(p == 0)=	-Inf;
x(p == 1)=	Inf;
end