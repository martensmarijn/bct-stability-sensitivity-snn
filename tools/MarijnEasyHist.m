function [bins, tarpf, fr] = MarijnEasyHist(cfg, firings, icell, exclneurons)

bins      = cfg.dtsc:cfg.dtsc:cfg.dtsc * cfg.Nt;
ii        = firings(:,2) > cfg.cci(icell) & firings(:,2) <= cfg.cci(icell + 1) & ~ismember(firings(:,2),exclneurons);
indxts    = round(firings(ii,1) / cfg.dtsc);
convt     = -3 * cfg.plot.sig_in_ms2:1:3 * cfg.plot.sig_in_ms2;
convy     = exp(-(convt .^ 2) ./ (2 * cfg.plot.sig_in_ms2 ^ 2))./sqrt(2 * pi * cfg.plot.sig_in_ms2 ^ 2);
convy     = convy / max(convy);
spikesbin = histc( indxts * cfg.dtsc, bins );
tarpf     = conv( spikesbin , convy ,'same');
fr        = sum(tarpf) / sum(convy) / (cfg.Nt * cfg.dtsc / 1000) / cfg.cci(icell + 1);