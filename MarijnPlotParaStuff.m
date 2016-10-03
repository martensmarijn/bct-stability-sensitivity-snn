function MarijnPlotParaStuff(cfg, inettyp, inet, ipara, MStore, RStore, SStore, optidx)

switch cfg.plot.options{optidx}
    case 'Figure3part1'
        if inet == 1 && inettyp == 3
            exclneurons = MStore(inettyp,inet,cfg.netpara).StimNeurons;
            figure(cfg.plot.icf(optidx)); hold on
            iparas = round(cfg.npara * [.22 .41 .6]);
            xlims = [find(diff(cfg.Ipatdrive)==1) * cfg.dtsc; find(diff(cfg.Ipatdrive)==1) * cfg.dtsc + 60];
            if ismember(ipara, iparas)
                firings     = RStore(inettyp,inet,ipara).stimfirings;
                jj          = find(firings(:,2) > cfg.cci(1) &  firings(:,2) <= cfg.cci(2)); 
                kk          = find(firings(:,2) > cfg.cci(2) &  firings(:,2) <= cfg.cci(3));
                ii          = find(firings(:,2) > cfg.cci(3) &  firings(:,2) <= cfg.cci(4));
                ss          = find(ismember(firings(:,2),exclneurons));
                subplot(2,2,find(ismember(iparas , ipara))); hold on
                    plot(firings(jj,1), firings(jj,2)                          ,'.','color',cfg.plot.nrnclr{1})
                    plot(firings(ss,1), firings(ss,2)                          ,'*','color',cfg.plot.nrnclr{1},'markersize',15)
                    plot(firings(kk,1), firings(kk,2) - cfg.cci(2) + cfg.nci(3),'.','color',cfg.plot.nrnclr{2})
                    plot(firings(ii,1), firings(ii,2) - cfg.cci(3)             ,'.','color',cfg.plot.nrnclr{3})
                    set(gca,cfg.plot.hp)
                    ylim([0 cfg.Ne])
                    xlim(xlims)
                    xlabel('time (ms)')
                    ylabel('neuron index')
            end
            dims   = size(squeeze(SStore(cfg.inettypes,1:cfg.Nnet,1:cfg.npara)));
            AmplAr = reshape([SStore.stimamp], dims);
            if mod(ipara,2) == 1
                subplot(2,2,4); hold on;
                    icell     = 1; % exc only
                    exclneurons =  MStore(inettyp,inet,cfg.netpara).StimNeurons;
                    firings = RStore(inettyp, inet, ipara).stimfirings;            
                    [bins, tarpfexclstim, ~] = MarijnEasyHist(cfg, firings, icell, exclneurons);
                    tarparar(:, icell, inettyp, inet, ipara) = tarpfexclstim;
                    spikedens = squeeze(tarparar(:,icell,inettyp,inet,ipara));
                    if find(squeeze(AmplAr(inettyp,inet,ipara)) < cfg.PeakThr)
                        plot(bins, spikedens,'m');
                    else
                        plot(bins, spikedens,'k');
                    end
                    set(gca,cfg.plot.hp); 
                    xlim(xlims)
                    xlabel('t (ms)');
                    ylabel('firing rate (Hz)');
            end
        end
    case 'Figure1part4'
        if inet == 1 && inettyp == cfg.inettypes(1) && ipara == 1
            a = [.02;  .1; .02];
            b = [ .2;  .2; .25]; 
            c = [-65; -65; -65];
            deltsyn     = [1 / cfg.tAMPA; 1 / cfg.tGABA; 1 / cfg.tGABA];
            expdeltsyn  = exp( -deltsyn * cfg.dtsc);
            Erev        = [cfg.EAMPA; cfg.EGABA; cfg.EGABA];
            Isyn        = 0;
            Ibg         = [2.45 3.5 .3];
            idx         = 0;
            timest      = [1000 1100] ./ cfg.dtsc;
            time        = 0:cfg.dtsc:(diff(timest) * cfg.dtsc);
            for posttype = 1:3
                for pretype = 1:3
                    idx = idx + 1;
                    v = c(posttype);
                    u = -13;
                    con         = cfg.SynConMat(posttype,pretype);
                    for t = 1:timest(end)
                        Isyn    = expdeltsyn(pretype) .* Isyn;
                        if t == timest(1) + 20 / cfg.dtsc 
                            Isyn = 1;
                        end
                        dum     = min(v, 50);
                        I       = Ibg(posttype) + con * (Erev(pretype) - dum) * Isyn;
                        v = v + cfg.dtsc * 0.5 * (0.04 * v .^2 + 5 * v + 140 - u + I);  % step 0.5
                        v = v + cfg.dtsc * 0.5 * (0.04 * v .^2 + 5 * v + 140 - u + I);  % for numerical stability
                        u = u + cfg.dtsc * a(posttype) .* (b(posttype) .* v - u);
                        V(t) = v;
                    end
                    subplot(3,3,idx); hold on
                    plot(time , V(timest(1):timest(2)))
                    ylim([-67 -65])
                    if posttype == 3
                        ylim([-64 -63])
                    elseif posttype == 2
                        ylim([-64.5 -62.5])
                    end
                end
            end
        end
end
