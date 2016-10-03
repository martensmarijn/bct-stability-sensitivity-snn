function MarijnPlotMetaStuff(cfg, MStore, SStore, RStore, optidx)

switch cfg.plot.options{optidx}
    case 'Figure1part2'
        figure(cfg.plot.icf(optidx)); hold on
        MeanDAr	= reshape([MStore.MeanD],size(MStore));
        clrs   	= ones(3,sqrt(length(cfg.var.Ne))) .* [ones(3,1) * linspace(0,.8,sqrt(length(cfg.var.Ne)))];
        subplot(3,2,1); hold on
            inettyp = 1;
            Nes = unique(cfg.var.Ne);
            MeanD = squeeze(MeanDAr(inettyp,:,:));
            for i = 1:length(Nes)
                idx = cfg.var.Ne == Nes(i);
                errorbar(cfg.var.pconSyn(idx), mean(MeanD(:,idx)), std(MeanD(:,idx)) ./ sqrt(length(MeanD(:,idx))), 'color',clrs(:,i))
                lgndstr{i} = ['N = ' num2str(Nes(i))]; %#ok<*AGROW>
            end
            legend(lgndstr)
            set(gca, cfg.plot.hp);
         subplot(3,2,2); hold on
            inettyp = 1;
            pcons = unique(cfg.var.pconSyn);
            MeanD = squeeze(MeanDAr(inettyp,:,:));
            for i = 1:length(pcons)
                idx = cfg.var.pconSyn == pcons(i);
                errorbar(cfg.var.Ne(idx), mean(MeanD(:,idx)), std(MeanD(:,idx)) ./ sqrt(length(MeanD(:,idx))), 'color',clrs(:,i))
                lgndstr{i} = ['pcon = ' num2str(pcons(i))];
            end
            legend(lgndstr)
            set(gca, cfg.plot.hp);
         subplot(3,2,3); hold on
            for inettyp = cfg.inettypes
                MeanD = squeeze(MeanDAr(inettyp,:,:));
                idx = abs(cfg.var.pconSyn - 0.05) < 0.01;
                errorbar(cfg.var.Ne(idx), mean(MeanD(:,idx)), std(MeanD(:,idx)) ./ sqrt(length(MeanD(:,idx))), 'color',cfg.plot.clr{inettyp})
                meandtt(inettyp,:,:) = MeanD(:,idx);
                lgndstr{inettyp} = num2str(inettyp);
            end
            plot(cfg.var.Ne(1:5),ttest2(squeeze(meandtt(1,:,:)),squeeze(meandtt(2,:,:))))
            ptitle = cfg.var.pconSyn(idx);
            title(['pcon = ' num2str(ptitle(1))])
            legend(cfg.plot.titlar(cfg.inettypes))
            set(gca, cfg.plot.hp);
         for ii = 1:3
             pvals = [0.01 0.05 0.15];
             subplot(3,2,3+ii); hold on
                inettyp = 1;
                MeanD = squeeze(MeanDAr(inettyp,:,:));
                idx = abs(cfg.var.pconSyn - pvals(ii)) < 0.01; %matlab precision
                MeanDPCOR = mean(MeanD(:,idx));
                for inettyp = cfg.inettypes
                    MeanD = squeeze(MeanDAr(inettyp,:,:));
                    errorbar(cfg.var.Ne(idx), mean(MeanD(:,idx)) ./ MeanDPCOR, std(MeanD(:,idx)) ./ sqrt(length(MeanD(:,idx))), 'color',cfg.plot.clr{inettyp})
                end
                ptitle = cfg.var.pconSyn(idx);
                title(['pcon = ' num2str(ptitle(1))])
                legend(cfg.plot.titlar(cfg.inettypes))
                set(gca, cfg.plot.hp);
         end  
	case 'Figure1part3'
        figure(cfg.plot.icf(optidx)); hold on
        DEGS = []; IDS = []; ODS = []; clrorder = [];
        for inettyp = cfg.inettypes
            ID{inettyp} = []; OD{inettyp} = []; DEG{inettyp} = [];
            for inet = 1:cfg.Nnet
               [id,od,deg] = degrees_dir(full(MStore(inettyp,inet,1).M(1:cfg.Ne,1:cfg.Ne)));
               ID{inettyp}   = [ID{inettyp} id];
               OD{inettyp}   = [OD{inettyp} od];
               DEG{inettyp}  = [DEG{inettyp} deg];
            end
            DEGS = [DEGS; DEG{inettyp}];
            IDS  = [IDS; ID{inettyp}];
            ODS  = [ODS; OD{inettyp}];
            clrorder = [clrorder; cfg.plot.clrrgb(inettyp,:)];
        end
        subplot(3,1,1); hold on; title('sum degrees'); set(gca,'colororder',clrorder);
            nbars = max(max(DEGS)) - min(min(DEGS));
            [N,I] = hist(DEGS',nbars);
            sumn = mean(sum(N,1));
            plot(I,N ./ sumn); 
            legend(cfg.plot.titlar(cfg.inettypes))
        subplot(3,1,2); hold on; title('in degrees'); set(gca,'colororder',clrorder);
            nbars = max(max(IDS)) - min(min(ODS));
            [N,I] = hist(IDS',nbars);
            plot(I,N ./ sumn)
        subplot(3,1,3); hold on; title('out degrees'); set(gca,'colororder',clrorder);
            nbars = max(max(ODS)) - min(min(ODS));
            [N,I] = hist(ODS',nbars);
            plot(I,N ./ sumn)
    case 'Figure2part2'
        figure(cfg.plot.icf(optidx)); hold on
            dimsS   = [length(cfg.inettypes) cfg.Nnet cfg.npara];
            dimsM   = [length(cfg.inettypes) cfg.Nnet];
            frA     = reshape([SStore(cfg.inettypes,1:cfg.Nnet,1:cfg.npara).fr], dimsS);            
            MeanDA	= reshape([MStore(cfg.inettypes,1:cfg.Nnet).MeanD], dimsM);
            id = 0;
            for inettyp = cfg.inettypes
                id = id + 1; fr = []; MeanD = []; kn = [];
                for inet = 1:cfg.Nnet
                    knct(inet,:) = MStore(inettyp,inet,1).kncentr;
                    for ipara = 1:cfg.npara %11
                        fr       = [fr frA(inettyp,inet,ipara)];
                        MeanD    = [MeanD MeanDA(inettyp,inet)];
                        kn       = [kn find(diff(MStore(inettyp,inet,1).kncentr < 1))];
                    end
                end
                subplot(3,2,1); hold on; xlabel('kcore (number of nodes in core)'); ylabel('fraction of neurons in k-core (%)')
                    krange = 20:45;
                    KN = squeeze(knct(:,krange));
                    KN = KN / max(max(KN)) * 100;
                    KNtt{inettyp} = KN;
                    errorbar(1:size(KN,2), mean(KN,1), std(KN,1) / sqrt(size(KN,1)), 'color',cfg.plot.clr{inettyp})
                    set(gca,'xtick',1:length(krange),'xticklabel',krange)
                [Xfrsort, Idx]  = sort(fr);
                YMeanDsort    	= MeanD(Idx);
                Yknsort         = kn(Idx);
                [pMeanD, errMeanD] = polyfit(Xfrsort, YMeanDsort, 1);
                [pkn, errkn]   	   = polyfit(Xfrsort, Yknsort, 1);
                [~, PcorrMeanD] = corrcoef(Xfrsort, YMeanDsort);
                [~, Pcorrkn]    = corrcoef(Xfrsort, Yknsort);
                y_fitMeanD     	= polyval(pMeanD, Xfrsort, errMeanD);
                y_fitkn         = polyval(pkn, Xfrsort, errkn);
                y_difMeanD    	= YMeanDsort - y_fitMeanD;
                y_difkn         = Yknsort - y_fitkn;
                SSdifMeanD    	= sum(y_difMeanD .^ 2);
                SSdifkn         = sum(y_difkn .^ 2);
                SStotMeanD  	= (length(YMeanDsort) - 1) * var(YMeanDsort);
                SStotkn         = (length(Yknsort) - 1) * var(Yknsort);
                R_MeanD         = 1- SSdifMeanD / SStotMeanD;
                R_kn            = 1  - SSdifkn / SStotkn;
                subplot(3,2,2); hold on; ylabel('mean pair distance'); xlabel('firing rate (Hz)')
                    h_MeanD(id) = plot(Xfrsort, y_fitMeanD,'color',cfg.plot.clr{inettyp},'linewidth',3);  
                    plot(Xfrsort, YMeanDsort, '.','color',cfg.plot.clr{inettyp})
                    lgndstr_MeanD{id} = ['R: ' num2str(R_MeanD) ' P: ' num2str(PcorrMeanD(2,1))];
                subplot(3,2,3); hold on; ylabel('neurons in k-core'); xlabel('firing rate (Hz)')
                    h_kn(id)   = plot(Xfrsort, y_fitkn, 'color',cfg.plot.clr{inettyp},'linewidth',3);  
                    plot(Xfrsort, Yknsort, '.','color',cfg.plot.clr{inettyp})
                    lgndstr_kn{id} = ['R: ' num2str(R_kn) ' P: ' num2str(Pcorrkn(2,1))];
            end
            legend(h_MeanD,lgndstr_MeanD)
            legend(h_kn,lgndstr_kn)
            H = ttest2(KNtt{1}(:,:),KNtt{2}(:,:));
            subplot(3,2,1); plot(H);
   	case 'Figure2part3'
        figure(cfg.plot.icf(optidx)); hold on
       
            dimsS     = [length(cfg.inettypes) cfg.Nnet cfg.npara];
            burstrate = reshape([SStore(cfg.inettypes,1:cfg.Nnet,1:cfg.npara).burstrate], dimsS);
            frA       = reshape([SStore(cfg.inettypes,1:cfg.Nnet,1:cfg.npara).fr], dimsS);
            if cfg.Nnet > 1 && cfg.npara > 1
                for inettyp = cfg.inettypes
                    if length(cfg.inettypes) == 1
                        fr      = squeeze(frA);
                        brate   = squeeze(burstrate);
                    else
                        fr      = squeeze(frA(inettyp,:,:));
                        brate   = squeeze(burstrate(inettyp,:,:));
                    end
                    subplot(2,1,1); hold on; xlabel('noise'); ylabel('mean firing rate')
                        errorbar(cfg.var.noisehz, mean(fr,1), std(fr,1) / sqrt(size(brate,1)), 'color',cfg.plot.clr{inettyp})
                        plot(cfg.var.noisehz, mean(fr,1), 'color', cfg.plot.clr{inettyp})
                        frtt{inettyp} = fr;
                        if inettyp == cfg.inettypes(end); legend(cfg.plot.titlar(cfg.inettypes)); end
                    subplot(2,1,2); hold on; xlabel('noise'); ylabel('burst rate')
                        errorbar(cfg.var.noisehz, mean(brate,1),std(brate,1) / sqrt(size(brate,1)), 'color', cfg.plot.clr{inettyp})
                        plot(cfg.var.noisehz, mean(brate,1), 'color', cfg.plot.clr{inettyp})
                        bratett{inettyp} = brate;
                end
            end
            H = ttest2(frtt{1},frtt{2});
            subplot(2,1,1); plot(cfg.var.noisehz,H)
            H = ttest2(bratett{1},bratett{2});
            subplot(2,1,2); plot(cfg.var.noisehz,H)
    case 'Figure2part4'
        figure(cfg.plot.icf(optidx)); hold on
        inettypes = [1 2];
        for inettyp = inettypes
            frse_all = []; InD_all = []; OutD_all = []; InW_all = [];
            ipara = cfg.npara;
            cfg = MarijnBuildFinalCfg(cfg, ipara); %get cfg.netpara right
            for inet = 1:cfg.Nnet
                con  = MarijnMakeCon(cfg,inettyp,inet,MStore);
                con  = con(1:cfg.Ne,1:cfg.Ne);
                [Pf(inet), Pw(inet), ~, ~, ~] = MarijnCalcCorrelationStr(con);
                frse_all = [frse_all; SStore(inettyp, inet, ipara).frse];
                OutD_all = [OutD_all; sum(con > 0, 1)'];
                InD_all  = [InD_all; sum(con > 0, 2)];
                InW_all  = [InW_all; sum(con, 2)];
                if inet == 1 && inettyp == 2
                    Econ = nan(cfg.Ne); 
                    Econ(con > 0) = con(con > 0);
                    subplot(length(inettypes) + 2, 2 , inettypes(end) * 2 + 1); hold on; xlabel('w'); ylabel('p');
                        edges = linspace(0,cfg.SynConMat(1,1) * 59,60);
                        N = histc(nansum(Econ), edges);
                        bh = bar(edges,N);
                        set(bh,'edgecolor','none');
                        xlim([-edges(2) edges(end) + edges(2)])
                    subplot(length(inettypes) + 2, 2 , inettypes(end) * 2 + 3); hold on; xlabel('w'); ylabel('p');
                        edges = linspace(0,0.02,40);
                        N = histc(Econ(:), edges);
                        bh = bar(edges,N);
                        set(bh,'edgecolor','none');
                        xlim([-edges(2) edges(end) + edges(2)])
                end
            end
            frse_all = frse_all + rand(length(frse_all),1);
            subplot(length(inettypes) + 2,2,inettyp * 2 - 1); hold on; ylabel('firing rate (Hz)'); xlabel('in-degree')
                [XD_all, Ix_all] = sort(InD_all);
                Y_all       = frse_all(Ix_all);
                [p_all,err_all] = polyfit(XD_all,Y_all,1);
                y_fitD_all   = polyval(p_all,XD_all,err_all);
                plot(InD_all(1:480),frse_all(1:480),'.k')
                plot(XD_all,y_fitD_all,'k');
            subplot(length(inettypes) + 2,2,inettyp * 2); hold on; ylabel('firing rate (Hz)'); xlabel('out-degree')
                [XD_all, Ix_all] = sort(OutD_all);
                Y_all       = frse_all(Ix_all);
                [p_all,err_all] = polyfit(XD_all,Y_all,1);
                y_fitD_all   = polyval(p_all,XD_all,err_all);
                plot(OutD_all(1:480),frse_all(1:480),'.k')
                plot(XD_all,y_fitD_all,'k');
            if inettyp == 2
                subplot(length(inettypes) + 2,2,inettyp * 2 + 2); hold on; ylabel('degree of ACOR windows');
                    edges = linspace(-1,-.2,41);
                    N = histc(Pw, edges);
                    bh(1) = bar(edges,N);
                    N = histc(Pf, edges);
                    bh(2) = bar(edges,N);
                    xlim([-1.05 -0.15])
                    set(bh(1),'facecolor','k','edgecolor','none','barwidth',1);
                    set(bh(2),'facecolor','r','edgecolor','none','barwidth',.5);
                    legend(bh,'PCOR/ACOR - 1','slope')
            end
        end
        inettyp = 1;
        ipara = 4; %12
        frA       = reshape([SStore(inettyp,1:cfg.Nnet,ipara).fr], [1 cfg.Nnet 1]);
        for inet = 1:cfg.Nnet
            con  = MarijnMakeCon(cfg,inettyp,inet,MStore);
            con  = con(1:cfg.Ne,1:cfg.Ne);
            [Pf(inet), Pw(inet), ~, ~, ~] = MarijnCalcCorrelationStr(con);
        end
%         subplot(4,2,7); hold on; ylabel('slope'); xlabel('fr');
%             plot(frA,Pf,'.')
        subplot(4,2,8); hold on; ylabel('PCOR/ACOR - 1'); xlabel('fr');
            plot(frA,Pw,'.')
    case 'Figure4part3'
        figure(cfg.plot.icf(optidx)); hold on
            if isfield(cfg.var,'NFLEX_stim')
                vars = cfg.var.NFLEX_stim;
                titler = 'Nstim';
            elseif isfield(cfg.var,'SynConMat')
                vars = cfg.var.SynConMat;
                titler = 'Wee';
            elseif isfield(cfg.var,'noisehz')
                vars = cfg.var.noisehz;
                titler = 'Ibg (Hz)';
            end
            dims = [length(cfg.PeakThrROC) length(cfg.inettypes) cfg.Nnet cfg.npara];
            TP   = reshape([SStore(cfg.inettypes,1:cfg.Nnet,1:cfg.npara).tp], dims);
            FP   = reshape([SStore(cfg.inettypes,1:cfg.Nnet,1:cfg.npara).fp], dims);
            for inettyp = 1:length(cfg.inettypes)
                clr    = cfg.plot.clr{cfg.inettypes(inettyp)};           
                for ipara = 1:cfg.npara
                    for inet = 1:cfg.Nnet
                        tp = squeeze(TP(:,inettyp,inet,ipara));
                        fp = squeeze(FP(:,inettyp,inet,ipara));
                        tpexample(inet,:) = tp;
                        fpexample(inet,:) = fp;
                        x      = flipud([1; fp]);
                        y      = flipud([1; tp]);
                        xdiff  = [x(1); diff(x)];
                        auc1   = sum(y .* xdiff);
                        auc2   = sum([0; y(1:end - 1)] .* xdiff);
                        AUCnets(inet) = mean([auc1, auc2]);
                    end
                    AUC(:,ipara) = AUCnets;
                    plotparas = [1 7 11];
                    if ismember(ipara,plotparas)
                        subplot(2,2,find(plotparas == ipara)); hold on
                            plot(mean(fpexample,1),mean(tpexample,1),'color',clr,'linewidth',2)
                            title([titler ' = ' num2str(vars(ipara))])
                            xlabel('false positive rate')
                            ylabel('true positive rate')
                    end
                end
                subplot(2,2,4); hold on; ylim([.45 1])
                    statAUC{inettyp} = AUC;
                    errorbar(vars, nanmean(AUC), nanstd(AUC) / sqrt(size(AUC,1)),'color',clr)
                    xlabel(titler)
                    ylabel('AUC')
            end
            [H P] = ttest2(statAUC{1},statAUC{2});
            plot(vars,H,'k')
            subplot(2,2,1)
            legend(cfg.plot.titlar(cfg.inettypes),'location','southeast')
    case 'Figure6part1'
        figure(cfg.plot.icf(optidx)); hold on
        fname   = fieldnames(cfg.var);
        xvar    = getfield(cfg.var,fname{1});
        inettyp = 3;
        inets   = 2%[4 5]% 43:45;
        for ipara = 1:cfg.npara;
            cfg = MarijnBuildFinalCfg(cfg, ipara); %get cfg.netpara right
            for inet = 1:cfg.Nnet
                con = RStore(inettyp,inet,ipara).constdp(1:cfg.Ne,1:cfg.Ne);
                conall = con;
                [~,Istdp] = sort(con(:));
                con(Istdp(1:round(0.95 * cfg.Ne^2))) = 0;
                [Pf(ipara,inet), Pw(ipara,inet), ~, ~, ~] = MarijnCalcCorrelationStr(con);
                frse = SStore(inettyp, inet, ipara).frse;
                InD = sum(con > 0, 2);
                    Y       = frse;
                    [XD, Ix] = sort(InD);
                    Y       = Y(Ix);
                    p = polyfit(XD,Y,1);
                    slopeD(ipara,inet) = p(1);
                InW = sum(con, 2);
                    Y       = frse;
                    [XW, Ix] = sort(InW);
                    Y       = Y(Ix);
                    p = polyfit(XW,Y,1);
                    slopeW(ipara,inet) = p(1);
                if ipara == cfg.npara
                    inm = ismember(inets,inet);
                    if sum(inm) == 1
                        subplot(2,4,1 + find(inets == inet)); hold on; title(['after STDP']); xlabel('in degree'); ylabel('out degree');
                            [~, ~, y_fit, sortInD, sortOutD] = MarijnCalcCorrelationStr(con);
                            plot(sortInD, y_fit, 'color',cfg.plot.clr{inettyp},'linewidth',2);
                            plot(sortInD, sortOutD,'.','color',cfg.plot.clr{inettyp}) 
                            axis([0 60 0 60])
                            axis square
                            M = MStore(inettyp,ipara).M(1:cfg.Ne,1:cfg.Ne);
                            [~, ~, y_fit, sortInD, sortOutD] = MarijnCalcCorrelationStr(M);
                            subplot(2,4,1); hold on; title('UCOR'); xlabel('in degree'); ylabel('out degree'); axis([0 100 0 100]); axis square
                                plot(sortInD, y_fit, 'color',cfg.plot.clr{inettyp},'linewidth',2);
                                plot(sortInD, sortOutD, '.','color',cfg.plot.clr{inettyp})
                            [~,Istpdall] = sort(conall(:));
                            remaining = conall(Istpdall(round(0.95 * cfg.Ne ^ 2) + 1:end));
                            pruned = conall(Istpdall(cfg.Ne ^ 2 - length(find(M)):round(0.95 * cfg.Ne ^ 2)));
                        if inm(end) == 1
                            subplot(2,4,4); hold on; xlabel('mean w in'); ylabel('p');
                                edges = linspace(0,cfg.SynConMat(1,1) * 59,60);
                                N = histc(nansum(con), edges);
                                bh = bar(edges,N,'g');
                            axis square
                                set(bh,'edgecolor','none');
                                xlim([-edges(2) edges(end) + edges(2)])                        	
			 end
                         subplot(2,4,4 + find(inets == inet)); hold on; xlabel('w'); ylabel('p');
                            edges = linspace(0,0.014,30);
                            Nremain = histc(remaining, edges);
                            Nprun = histc(pruned, edges);
                            bh(1) = bar(edges,Nremain,'g');
                            bh(2) = bar(edges,Nprun,'r');
                            set(bh(1),'edgecolor','none');
                            set(bh(2),'edgecolor','none');
                            axis square
                            % xlim([-edges(2) edges(end) + edges(2)])
                            % set(gca,'xtick',0:0.004:0.014)
                    end
                end
            end
        end
        xax = repmat(xvar',1,cfg.Nnet);
        subplot(2,4,7); hold on; ylabel('degree of ACOR windows'); xlabel(fname)
            edges = linspace(-1.2,-.2,41);
            N = histc(Pw, edges);
            bh(1) = bar(edges,N);
            N = histc(Pf, edges);
            bh(2) = bar(edges,N);
            set(bh(1),'facecolor','k','edgecolor','none','barwidth',1);
            set(bh(2),'facecolor','r','edgecolor','none','barwidth',.5);
            legend(bh,'PCOR/ACOR - 1','slope')
                            axis square
            xlim([-1.20 -0.3])
      %  subplot(3,3,8); hold on; ylabel('slope indegree vs frs'); xlabel(fname);
      %      plot(xax, slopeD,'k');
      %  subplot(3,3,9); hold on; ylabel('slope indegree vs mean presyn str'); xlabel(fname)
       %     plot(xax, slopeW,'b');
end
