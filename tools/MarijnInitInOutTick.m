function [InTick, OutTick] = MarijnInitInOutTick(inettyp, N, mu, th, sigx, sigy)

sigscal     = mean(mu) / 3;
R           = [cos(th) sin(th); -sin(th) cos(th)];
InDeg       = 3 * (N + abs(diff(mu * N)));
OutDeg      = 0;

if ismember(inettyp,1:3)
    while abs(sum(InDeg) - sum(OutDeg)) > abs(diff(mu * N))
        x       = R * diag([sigx sigy]) * randn(2, N * 10);
        x       = round(repmat([mu(1); mu(2)],[1 N * 10]) + sigscal * x);
        ii      = find(x(1,:) >= 1 & x(1,:) <= 2 * mu(1) & x(2,:) >= 1 & x(2,:) <= 2 * mu(2));
        jj      = randperm(length(ii));
        ii      = ii(jj(1:N));
        InDeg   = x(1,ii);
        OutDeg  = x(2,ii);
    end
    if inettyp == 3
        OutDeg = OutDeg(randperm(length(OutDeg)));
    end
end
if ismember(inettyp,4)
    while abs(sum(InDeg) - sum(OutDeg)) >  abs(diff(mu * N))
        R1      = [cos(th) sin(th); -sin(th) cos(th)];
        R2      = [cos(-th) sin(-th); -sin(-th) cos(-th)];
        x       = [R1 * diag([sigx sigy]) * randn(2, N * 5) R2 * diag([sigx sigy]) * randn(2, N * 5)];
        x       = round(repmat([mu(1); mu(2)],[1 N * 10]) + sigscal * x);
        ii      = find(x(1,:) >= 1 & x(1,:) <= 2 * mu(1) & x(2,:) >= 1 & x(2,:) <= 2 * mu(2));
        jj      = randperm(length(ii));
        ii      = ii(jj(1:N));
        InDeg   = x(1,ii);
        OutDeg  = x(2,ii);
    end
end

Excess = sum(OutDeg) - sum(InDeg) - diff(mu * N); % we add the lowest and substract from the highest
%to get the sums the same, which is a requirement
if Excess > 0
    dum1                = round(Excess / 2);
    dum2                = Excess - dum1;
    [Nothing Indx]      = sort(OutDeg,'Descend');
    OutDeg(Indx(1:dum1))= OutDeg(Indx(1:dum1)) - 1;
    [Nothing Indx]      = sort(InDeg,'Ascend');
    InDeg(Indx(1:dum2)) = InDeg(Indx(1:dum2)) + 1;    
else
    dum1                = round( abs(Excess) / 2);
    dum2                = abs(Excess) - dum1;
    [Nothing Indx]      = sort(OutDeg,'Ascend');
    OutDeg(Indx(1:dum1))= OutDeg(Indx(1:dum1)) + 1;
    [Nothing Indx]      = sort(InDeg,'Descend');
    InDeg(Indx(1:dum2)) = InDeg(Indx(1:dum2)) - 1;
end

dum          = 1 + [0 cumsum(InDeg(1:end - 1))]; % Vector of indegree sums padded by one minus the last
InTick       = zeros(sum(InDeg), 1); % Vector of incoming edges
InTick(dum)  = 1;
InTick       = cumsum(InTick); % Cumulative

dum     	 = 1 + [0 cumsum(OutDeg(1:end - 1))];
OutTick      = zeros(sum(OutDeg), 1);
OutTick(dum) = 1;
OutTick      = cumsum(OutTick);
OutTick      = OutTick(randperm(length(OutTick))); %randomize