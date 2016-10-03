function cfg = MarijnBuildFinalCfg(cfg, ipara)

fnames = fieldnames(cfg.var);
for isel = 1:length(fnames)
    params = cfg.var.(fnames{isel});
    cfg.(fnames{isel})(cfg.iparidx(1), cfg.iparidx(2)) = params(ipara);
end

cfg.Ni1      = round(0.25 * cfg.Ne * 0.5); % 10% FS PV+ target inhibitory cells
cfg.Ni2      = round(0.25 * cfg.Ne * 0.5); % 10% LTS SOM+ target PC cells
cfg.Ntot     = cfg.Ne + cfg.Ni1 + cfg.Ni2;
cfg.nci      = [cfg.Ne cfg.Ni1 cfg.Ni2];
cfg.cci      = cumsum([0 cfg.nci]);

if cfg.stdpselect_t > cfg.Nt; cfg.stdpselect_t = cfg.Nt; end

if ipara == 1 || isfield(cfg.var,'Ne') || isfield(cfg.var,'Ni1') || isfield(cfg.var,'Ni2') || isfield(cfg.var,'pconSyn')
    cfg.netpara = ipara;
end