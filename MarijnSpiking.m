function RStore = MarijnSpiking(cfg, inettyp, inet, ipara, MStore, RStore, con)

rng(cfg.vary_spike_with_inet * inet + 100 * cfg.whattype, 'twister'); % vary seed with Net

a = [.02     * ones(cfg.Ne,1) ;  .1     * ones(cfg.Ni1,1) ; .02     * ones(cfg.Ni2,1) ];
b = [ .2     * ones(cfg.Ne,1) ;  .2     * ones(cfg.Ni1,1) ; .25     * ones(cfg.Ni2,1) ]; 
c = [-65 + 5 * randn(cfg.Ne,1); -65 + 5 * randn(cfg.Ni1,1); -65 + 5 * randn(cfg.Ni2,1)];
d = [  8     * ones(cfg.Ne,1) ;   2     * ones(cfg.Ni1,1) ;   2     * ones(cfg.Ni2,1) ];
v = c;
u = 3 * randn(cfg.Ntot,1);

Isyn       = zeros(cfg.Ntot,1);
deltsyn    = [ones(cfg.Ne,1) / cfg.tAMPA; ones(cfg.Ni1 + cfg.Ni2,1) / cfg.tGABA];
Erev       = [cfg.EAMPA * ones(cfg.Ne,1); cfg.EGABA * ones((cfg.Ni1 + cfg.Ni2),1)];
expdeltsyn = exp( - cfg.dtsc * deltsyn);
taustdp    = exp( - cfg.dtsc * 1 / cfg.taustdp );
STDPD      = sparse(zeros(cfg.Ntot));
STDPP      = sparse(zeros(cfg.Ntot));
w0         = cfg.SynConMat(1,1) / cfg.w0fac;
STDPl      = 0.1;
STDPtau    = 20; %in ms
STDPa      = 0.11;
D          = sparse(zeros(cfg.Ne));
P          = sparse(zeros(cfg.Ne));

if ~isfloat(con)
    con = MarijnMakeCon(cfg,inettyp,inet,MStore);
end
constdp         = [];
[ia ja ka]      = find(con);
repeat          = true;
needsbasefire   = cfg.needsbasefire;

while repeat
    rng(cfg.vary_spike_with_inet * inet + 100 * cfg.whattype, 'twister'); % reset seed for every repeat
    firings         = [];
    repeat          = false;
   	Ipatdrive       = cfg.Ipatdrive;
    statetimes      = [1; find(diff(Ipatdrive) == 1) - 2; -1];
    nistim          = MStore(inettyp, inet, cfg.netpara).StimNeurons;
    st              = 1;
    store_basefire  = false;
    if cfg.only1basefire
        xpara = 1;
    else
        xpara = ipara;
    end
    if needsbasefire && ipara == xpara
        needsbasefire  = false;
        Ipatdrive      = Ipatdrive .* 0;
        store_basefire = true; 
        if cfg.patdrive > 0
            repeat     = true;
        end
    end
    for t = 1:cfg.Nt % MAIN LOOP %
        %%%%%%%%%%%%%%%%%%%%%%
        %   state recall     %
        %%%%%%%%%%%%%%%%%%%%%%
        if cfg.dostdp
            [ia ja ka]  = find(con);
        end
        if t == statetimes(st)
            IpatNe = zeros(cfg.Ne, 1);
            if cfg.patdrive > 0
                IpatNe(nistim(randperm(length(nistim),cfg.NFLEX_stim))) = cfg.CurFlex;
                if store_basefire
                    RStore(inettyp,xpara).store{st}.v    = v;
                    RStore(inettyp,xpara).store{st}.u    = u;
                    RStore(inettyp,xpara).store{st}.Isyn = Isyn;
                else
                    if cfg.needsbasefire == true
                        v    = RStore(inettyp,xpara).store{st}.v;
                        u    = RStore(inettyp,xpara).store{st}.u;
                        Isyn = RStore(inettyp,xpara).store{st}.Isyn;
                    end
                end
                RStore(inettyp,ipara).store{st}.IpatNe = IpatNe;
            end
            st = st + 1;
        end
        %%%%%%%%%%%%%%%%%%%%%%
        %   input current    %
        %%%%%%%%%%%%%%%%%%%%%%
        IstimN = zeros(cfg.Ne, 1);
        IstimN(randi(cfg.Ne, round(2 * rand * cfg.noisehz * cfg.dtsc * cfg.Ne),1)) = cfg.IstimAmpRandom;
        I = [IpatNe  * Ipatdrive(t) + cfg.Ibg(1) + cfg.Isd(1) * randn(cfg.Ne, 1) + IstimN; 
                                      cfg.Ibg(2) + cfg.Isd(2) * randn(cfg.Ni1,1)         ;
                                      cfg.Ibg(3) + cfg.Isd(3) * randn(cfg.Ni2,1)         ];
        %%%%%%%%%%%%%%%%%%%%%%
        %   cell stuff       %
        %%%%%%%%%%%%%%%%%%%%%%
        fired   = find(v >= 30);
        firings = [firings; cfg.dtsc * t + 0 * fired, fired];
        v(fired)= c(fired);
        u(fired)= u(fired) + d(fired);
        Isyn    = expdeltsyn .* Isyn;
        Isyn(fired) = Isyn(fired) + 1;
        dum     = min(v, 50);
        I       = I + sparse(ia, ja, ka .* (Erev(ja) - dum(ia)), cfg.Ntot, cfg.Ntot) * Isyn;
        v       = v + cfg.dtsc * 0.5 * (0.04 * v .^2 + 5 * v + 140 - u + I);
        v       = v + cfg.dtsc * 0.5 * (0.04 * v .^2 + 5 * v + 140 - u + I);
        u       = u + cfg.dtsc * a .* (b .* v - u);
        %%%%%%%%%%%%%%%%%%%%%%
        % plasticity stuff   %
        %%%%%%%%%%%%%%%%%%%%%%
        if cfg.dostdp 
            if mod(t, 3) == 1 && t * cfg.dtsc > 500 && ~isempty(firings) % network scaling
                ii          = firings(:,2) <= cfg.Ne;
                indxts      = firings(ii,1) > t * cfg.dtsc - 500;
                frcur       = length(find(indxts)) / 0.5 / cfg.Ne;
                frtgt       = cfg.frtar;
                dt          = 3 / (1000 / cfg.dtsc);
                con(1:cfg.Ne,1:cfg.Ne) = con(1:cfg.Ne,1:cfg.Ne) + con(1:cfg.Ne,1:cfg.Ne) .* (frtgt - frcur) * dt ./ cfg.homtau;
            end
            if cfg.morrisonstdp % morrison 2007 STDP
                if ~isempty(fired) && t < cfg.stdpselect_t
                    firede          = zeros(1,cfg.Ne,1);
                    firede(fired(fired < cfg.Ne)) = 1;
                    prespikes 	= ones(cfg.Ne,1) * firede;
                    postspikes 	= firede' * ones(1,cfg.Ne);
                    con(1:cfg.Ne,1:cfg.Ne)  = con(1:cfg.Ne,1:cfg.Ne) + prespikes .* P;
                    con(1:cfg.Ne,1:cfg.Ne)  = con(1:cfg.Ne,1:cfg.Ne) + postspikes .* D;
                    conmin            = 0;
                    con(con < conmin) = conmin; %avoid neg power
                    D = D - STDPa * STDPl * postspikes .* con(1:cfg.Ne,1:cfg.Ne);
                    D = D - cfg.dtsc * D / STDPtau;
                    P = P + STDPl * w0 ^ 0.6 .* con(1:cfg.Ne,1:cfg.Ne) .^ 0.4 .* prespikes; 
                    P = P - cfg.dtsc * P / STDPtau;
                    conmax            = 0.013;
                    conmin            = 0;
                    con(con < conmin) = conmin;
                    con(con > conmax) = conmax;
                end
            else % STDP
                STDPD       = taustdp * STDPD;
                STDPP       = taustdp * STDPP;
                if ~isempty(fired) && t < cfg.stdpselect_t % STDP
                    firede           = fired(fired < cfg.Ne);
                    con(firede, :)   = con(firede,:) + cfg.Pamp * STDPP(firede,:);
                    con(:, firede)   = con(:,firede) - cfg.Damp * STDPD(:,firede);
                    STDPD(firede, 1:cfg.Ne) = 1;
                    STDPP(1:cfg.Ne, firede) = 1;
                    STDPD            = STDPD .* (con > 0);
                    STDPP            = STDPP .* (con > 0);
                    conmax           = 0.013;
                    conmin           = 0;
                    con(con < conmin) = conmin;
                    con(con > conmax) = conmax;
                end
            end	
 
            if t == cfg.stdpselect_t % make connectivity 5%
                constdp = con(1:cfg.Ne,1:cfg.Ne);
                contemp = constdp;
                [~,Istdp] = sort(contemp(:));
                contemp(Istdp(1:round(0.95 * cfg.Ne^2))) = 0;
                con(1:cfg.Ne,1:cfg.Ne) = contemp;
                [ia ja ka]  = find(con);
            end
        end
    end
    if store_basefire
        RStore(inettyp,ipara).firings = firings;
    else
        RStore(inettyp,ipara).stimfirings = firings;
    end
end
if cfg.dostdp
    RStore(inettyp,ipara).con = con;
    RStore(inettyp,ipara).constdp = constdp;
end
