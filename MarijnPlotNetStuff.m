function MarijnPlotNetStuff(cfg, inettyp, inet, MStore, RStore, SStore, optidx)

switch cfg.plot.options{optidx}
	case 'Figure1part1'
        if inet == 1
            figure(cfg.plot.icf(optidx)); hold on
            M = MStore(inettyp, inet, :).M(1:cfg.Ne,1:cfg.Ne);
            subplot(3,2,inettyp)
                plot(sum(M > 0, 1) + rand(1,cfg.Ne),sum(M > 0, 2) + rand(cfg.Ne,1),'k.')
                set(gca, cfg.plot.hp);
                xlim([0 50]); ylim([0 50]);            
        end
  	case 'Figure3part2'    	
        if inet == 1 && inettyp == 2
            figure(cfg.plot.icf(optidx)); hold on
            dims = size(squeeze(SStore(cfg.inettypes,1:cfg.Nnet,1:cfg.npara)));   
            AmplAr              = reshape([SStore.stimamp], dims);
            dim                 = find(size(AmplAr) == cfg.npara) - 1;
            AmplArMeanNet       = squeeze(nanmean(AmplAr, dim))';
            AmplArThr           = AmplAr > cfg.PeakThr;
            AmplExceedThrArMean = 100 * squeeze(nanmean(AmplArThr, dim))';
            FirstPeakAr         = reshape([SStore.stimlat], dims);
            FirstPeakArMean     = squeeze(nanmean(FirstPeakAr, dim));
            for ipara = 1:cfg.npara               
                H(1,ipara) = ttest2(FirstPeakAr(1,:,ipara),FirstPeakAr(2,:,ipara));
                H(2,ipara) = ttest2(AmplAr(1,:,ipara),AmplAr(2,:,ipara));
                H(3,ipara) = ttest2(AmplArThr(1,:,ipara),AmplArThr(2,:,ipara));
            end
            H = H * 100;
            for inettyp = cfg.inettypes            
                subplot(3,1,1); hold on
                    h = plot(cfg.var.SynConMat, FirstPeakArMean(inettyp,:),'color', cfg.plot.clr{inettyp});
%                     FPA = squeeze(FirstPeakAr(inettyp,:,:));
%                     nonnan = ~(sum(isnan(FPA))==size(FPA,1));
%                     FPA = FPA(:,nonnan);
%                     xvars = cfg.var.SynConMat(nonnan);
%                     stdshade(FPA,0.15,cfg.plot.clr{inettyp});
%                     if inettyp == 1; plot(H(1,nonnan)); end
                    set(gca,cfg.plot.hp)
                    set(h,cfg.plot.h)
                    plot(cfg.var.SynConMat,H(1,:))
                    xlabel('Wee')
                    ylabel('peak latency (ms)')
                    xlim([.016 .04])
                subplot(3,1,2); hold on
                    h = plot(cfg.var.SynConMat, AmplArMeanNet(:,inettyp), 'color', cfg.plot.clr{inettyp});
                    set(gca,cfg.plot.hp)
                    set(h,cfg.plot.h)
                    plot(cfg.var.SynConMat,H(2,:))
                    xlabel('Wee')
                    ylabel('peak amplitude (Hz)')
                    xlim([.016 .04])
                subplot(3,1,3); hold on
                    h = plot(cfg.var.SynConMat, AmplExceedThrArMean(:,inettyp), 'color', cfg.plot.clr{inettyp});
                    set(gca,cfg.plot.hp)
                    set(h,cfg.plot.h)
                    plot(cfg.var.SynConMat,H(3,:))
                    xlabel('Wee')
                    ylabel('peak > threshold (%)')
                    xlim([.016 .04])
            end
            legend(cfg.plot.titlar(cfg.inettypes(cfg.inettypes)),'location','SouthEast')
        end
    case 'Figure2part1'
        if inet == 1 && inettyp == 3
            figure(cfg.plot.icf(optidx)); hold on
            times = cfg.plot.times;
            for ipara = 1:cfg.npara;   
                exclneurons = MStore(inettyp,inet,cfg.netpara).StimNeurons;         
                for icell = 1:length(cfg.nci) 
                    firings = RStore(inettyp, inet, ipara).firings;
                    [bins, tarpfexclstim, ~] = MarijnEasyHist(cfg, firings, icell, exclneurons);
                    tarparar(:, icell, inettyp, inet, ipara) = tarpfexclstim;
                end
            end
            xlabs = bins(times);
            offset = 12.5;
            for ipara = 1:2:cfg.npara;            
                spikedens   = squeeze(tarparar(:,:,inettyp,inet,ipara));
                for icell = [3 1 2]
                    plot(xlabs, spikedens(times,icell) + ipara * offset,'linewidth',2,'color',cfg.plot.nrnclr{icell});
                end
            end
        end  
    case 'Figure4part1'
        if inet == 6 && inettyp == 2
            figure(cfg.plot.icf(optidx)); hold on
            figure;inet = 11%inet = 9; 
	inettyp = 2; hold on
            
            offset 	= 12;
            times   = cfg.plot.times;
            ipara = 1;
            icell = 1;
            exclneurons = MStore(inettyp,inet,cfg.netpara).StimNeurons;
            firings = RStore(inettyp, inet, ipara).firings;
            [bins, spikedens, ~] = MarijnEasyHist(cfg, firings, icell, exclneurons);
            xalabs = bins(times);
            xlabs  = xalabs - xalabs(1);
            plot(xlabs, spikedens(times),'linewidth',2,'color',cfg.plot.nrnclr{icell});
            Idrive  = cfg.Ipatdrive(times);
            Idrive  = Idrive(1:1/cfg.dtsc:end);
            for ipara = 1:cfg.npara;            
                firings = RStore(inettyp,inet, ipara).stimfirings;
                [~, spikedens, ~] = MarijnEasyHist(cfg, firings, icell, exclneurons);
                plot(xlabs, spikedens(times) + offset * ipara,'linewidth',2,'color',cfg.plot.nrnclr{icell});
                stimon = find(diff(Idrive) == 1);
                for id = 1:length(stimon)
                    plot([stimon(id) stimon(id)], [0 offset * ipara],'k--')
                end
            end
            plot(Idrive * offset / 3 - offset / 1.5,'k')
            xlim([xlabs(1) xlabs(end)])
            xlabel('time (ms)')
            ylabel('firing rate (Hz)')
        end
    case 'Figure4part2'
        if inet == 5
            figure(cfg.plot.icf(optidx)); hold on
            thridx  = 3:length(cfg.PeakThrROC);
            dims    = [length(cfg.PeakThrROC) length(cfg.inettypes) cfg.Nnet cfg.npara];
            TP      = reshape([SStore(cfg.inettypes,1:cfg.Nnet,1:cfg.npara).tp], dims);
            FP      = reshape([SStore(cfg.inettypes,1:cfg.Nnet,1:cfg.npara).fp], dims);
            iparas  = [1 4 7];
            clrs    = {[0 0 0],[.5 .5 .5]};
            xlabs = cfg.PeakThrROC(thridx);
            if inettyp < 3
            for idx = 1:length(iparas)
                ipara = iparas(idx);
                subplot(length(cfg.inettypes) - 1,length(iparas),idx + 3 * (inettyp - 1)); hold on
                    tp = squeeze(TP(:,inettyp,:,ipara));
                    fp = squeeze(FP(:,inettyp,:,ipara));
                    tpmean = mean(tp,2);
                    tpsem  = std(tp,0,2) ./ size(tp,2);                   
                    fpmean = mean(fp,2);
                    fpsem  = std(fp,0,2) ./ size(fp,2); 
                    h(1) = bar(thridx,tpmean(thridx),'facecolor',clrs{2},'barwidth', 1,'edgecolor','none');
                    errorbar(thridx,tpmean(thridx),tpsem(thridx),'linestyle','none','color',clrs{2});
                    h(2) = bar(thridx,fpmean(thridx),'facecolor',clrs{1},'barwidth',.5,'edgecolor','none');
                    errorbar(thridx,fpmean(thridx),fpsem(thridx),'linestyle','none','color',clrs{1});
                    xlabel('threshold (Hz)')
                    ylabel('detection rate')
                    xlim([thridx(1) - 1 thridx(end) + 1])
                    ylim([0 1])
                    set(gca,'xtick',thridx(1:4:end),'xticklabel',xlabs(1:4:end))
                    title(['Nstim ' num2str(cfg.var.NFLEX_stim(ipara)) ' - ' cfg.plot.titlar{inettyp}])
            end
            legend(h, 'true positive rate','false positive rate')
            end
        end
end
