for inettyp = 1:3
    for idx = 1:1000
    sigx  	= 0.3;
    sigy  	= 1.0;
    mu(1)  	= cfg.pconSyn(1,1) * cfg.Ne;
    mu(2)   = cfg.pconSyn(1,1) * cfg.Ne;
    th = atan(mu(1) / mu(2));
    if inettyp == 2; 
        th = -th;
    end

    [InTick, OutTick]     = MarijnInitInOutTick(inettyp, cfg.Ne, mu, th, sigx, sigy);

    [~,Iuni] = unique(InTick * N + OutTick);
    Imulti(inettyp,idx)  = length(setdiff((1:length(InTick)), Iuni));
    Iself(inettyp,idx) = length(find(InTick == OutTick)'); % self-connections and double connections
    end
end

Imulti = Imulti / (cfg.Ne^2 * cfg.pconSyn(1,1)) * 100;
Iself = Iself / (cfg.Ne^2 * cfg.pconSyn(1,1)) * 100 ;

mIm = mean(Imulti')
mIs = mean(Iself')

sIm = std(Imulti')
sIs = std(Iself')

[hm,pm] = ttest2(Imulti(1,:),Imulti(2,:))
[hs,ps] = ttest2(Iself(1,:),Iself(2,:))

