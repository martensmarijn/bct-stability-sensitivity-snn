function MStore = MarijnGenNetwork(cfg, inettyp, inet, ipara, MStore)

if ipara == cfg.netpara
rng(inettyp + 5 * inet + 100 * ipara + 10000 * cfg.whattype,'twister');

sigx  	= 0.3;
sigy  	= 1.0;
mu(1)  	= cfg.pconSyn(1,1) * cfg.Ne;
mu(2)   = cfg.pconSyn(1,1) * cfg.Ne;
th = atan(mu(1) / mu(2));
if inettyp == 2; 
    th = -th;
end

if inettyp <= 4
    [InTickE, OutTickE]     = MarijnInitInOutTick(inettyp, cfg.Ne, mu, th, sigx, sigy);
    OutTickE                = MarijnPermuteDoubleCon(cfg.Ne, InTickE, OutTickE);
else
    Ne = round(cfg.Ne / 2);
    [InTickE_P, OutTickE_P] = MarijnInitInOutTick(1, Ne, mu, th, sigx, sigy);
    OutTickE_P              = MarijnPermuteDoubleCon(cfg.Ne, InTickE_P, OutTickE_P);
    [InTickE_A, OutTickE_A] = MarijnInitInOutTick(2, cfg.Ne - Ne, mu, -th, sigx, sigy);
    OutTickE_A              = MarijnPermuteDoubleCon(cfg.Ne, InTickE_A, OutTickE_A);
    if inettyp == 5
        InTickE          	= [InTickE_P; InTickE_A + cfg.Ne - Ne];
        OutTickE          	= [OutTickE_P; OutTickE_A + cfg.Ne - Ne];    
    elseif inettyp == 6
      	InTickE          	= [InTickE_A; InTickE_P + cfg.Ne - Ne];
        OutTickE          	= [OutTickE_A; OutTickE_P + cfg.Ne - Ne];
%     elseif inettyp == 7
%         InTickE         	= [InTickE_A;                InTickE_P + cfg.Ne - Ne];
%         OutTickE            = [OutTickE_P + cfg.Ne - Ne; OutTickE_A             ];
    end
    for iperm = 1:round(cfg.Ne * cfg.pconSyn(1,1)) * 20 % permute per neuron
        idx1 = []; idx2 = [];
        while length(idx1) < 1 || length(idx2) < 1
            idx1 = find(OutTickE == randi(Ne));
            idx2 = find(OutTickE == (randi(Ne) + (cfg.Ne - Ne)));
        end
        idxs = [idx1(randi(length(idx1))) idx2(randi(length(idx2)))];
        OutTickE(idxs) = OutTickE(fliplr(idxs));
    end
    OutTickE                = MarijnPermuteDoubleCon(cfg.Ne, InTickE, OutTickE);
end

M        = sparse(InTickE, OutTickE, ones(length(InTickE),1), cfg.Ntot, cfg.Ntot);
M        = M'; %Itskov FORMULAS ARE FOR OPPOSITE CONVENTION.
M        = full(M);

for i = 1:length(cfg.nci)
    indx{i} = 1:cfg.nci(i);
end
R = rand(cfg.Ntot);
for i = 1:length(cfg.nci) %receiver
    for j = 1:length(cfg.nci) %sender
        if ~[i == 1 && j == 1]
            M(cfg.cci(i) + indx{i}, cfg.cci(j) + indx{j}) = cfg.pconSyn(i,j) > R(cfg.cci(i) + indx{i}, cfg.cci(j) + indx{j});
        end
    end
end
OutDegree 	= sum(M(1:cfg.Ne,1:cfg.Ne) > 0, 1);
InDegree  	= sum(M(1:cfg.Ne,1:cfg.Ne) > 0, 2);
mn          = cfg.Ne * cfg.pconSyn(1,1);
Indx        = [];
rnge        = 0;
nstims      = max([10 cfg.NFLEX_stim]);
if isfield(cfg.var,'NFLEX_stim')
    nstims = max([cfg.var.NFLEX_stim nstims]);
end
while length(Indx) < nstims
    if inettyp <= 4
        Indx = find(OutDegree >= mn - rnge & OutDegree <= mn + rnge);
    else % stim the PCOR subset
        Indx = find(OutDegree(1:Ne) >= mn - rnge & OutDegree(1:Ne) <= mn + rnge);
    end
    rnge = rnge + 1;
end

MStore(inettyp,inet,cfg.netpara).M           = sparse(M);
MStore(inettyp,inet,cfg.netpara).OutDegree   = OutDegree;
MStore(inettyp,inet,cfg.netpara).InDegree    = InDegree;
MStore(inettyp,inet,cfg.netpara).StimNeurons = Indx(1:nstims);
end