function con = MarijnMakeCon(cfg,inettyp,inet,MStore)

M = MStore(inettyp,inet,cfg.netpara).M;
for i = 1:length(cfg.nci)
    indx{i} = 1:cfg.nci(i);
end
con = zeros(cfg.Ntot);
for i = 1:length(cfg.nci)
    for j = 1:length(cfg.nci)
        con(cfg.cci(i) + indx{i}, cfg.cci(j) + indx{j}) = cfg.SynConMat(i,j) * M(cfg.cci(i) + indx{i}, cfg.cci(j) + indx{j});
    end
end
if cfg.normweights
    scaling = (cfg.Ne * cfg.pconSyn(1,1) * cfg.SynConMat(1,1)) ./ sum(con(1:cfg.Ne,1:cfg.Ne),2) ;
    scalingmat = repmat(scaling,1,cfg.Ne);
    con(1:cfg.Ne,1:cfg.Ne) = con(1:cfg.Ne,1:cfg.Ne) .* scalingmat;
end