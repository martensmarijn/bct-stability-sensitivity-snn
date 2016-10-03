function OutTick = MarijnPermuteDoubleCon(N, InTick, OutTick)

last = 0;
Iabnorm = 1;
while ~isempty(Iabnorm)
    [~,Iuni] = unique(InTick * N + OutTick);
    Iabnorm  = unique([setdiff((1:length(InTick)), Iuni) find(InTick == OutTick)']); % self-connections and double connections

    if length(Iabnorm) ~= last % permute
        swap  = OutTick(Iabnorm);
        swap  = swap(randperm(length(swap)));
        OutTick(Iabnorm) = swap;
    else % loop stuck, permute like this
        Inorm = setdiff((1:length(InTick)), Iabnorm);
        Iswap = Inorm(randperm(length(Inorm),length(Iabnorm)));
        swap  = OutTick(Iswap);
        OutTick(Iswap)   = OutTick(Iabnorm);
        OutTick(Iabnorm) = swap;
    end
    last = length(Iabnorm);
end
