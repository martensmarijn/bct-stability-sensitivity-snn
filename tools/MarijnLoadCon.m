function con = MarijnLoadCon(cfg, inettyp, inet, datdir)

loadwhattype = cfg.loadwhattype;
loadwhatipara = cfg.loadwhatipara;
if loadwhattype > 0
    load(sprintf('%s/MarijnSpiking_%d.mat',datdir,loadwhattype));
    con = RStore(inettyp,inet,loadwhatipara).con;
else
    con = false;
end