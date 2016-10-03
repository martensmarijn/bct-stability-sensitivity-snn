function SStore = MarijnExtractSpikeStats(cfg, RStore, MStore)

for inettyp = cfg.inettypes
    for inet = 1:size(RStore,2)
        for ipara = 1:size(RStore,3)
            clearvars -except cfg RStore MStore SStore inettyp inet ipara
            cfg = MarijnBuildFinalCfg(cfg, ipara);
            cfg = MarijnGenStimulus(cfg);
            icell = 1;
            exclneurons = MStore(inettyp,inet,cfg.netpara).StimNeurons;
            %%%%%%%%%%%%%%%%%%%%%%
            % baseline responses %
            %%%%%%%%%%%%%%%%%%%%%%
            baseamps    = NaN;
            range       = NaN;
            burstrate   = NaN;
            frexcitatory= NaN;
            frse        = NaN;
            if cfg.needsbasefire
                if cfg.only1basefire
                    firings = RStore(inettyp, inet, 1).firings;
                else 
                    firings = RStore(inettyp, inet, ipara).firings;
                end
                frse = zeros(cfg.Ne,1);
                for ineuron = 1:cfg.Ne
%                     if ~ismember(ineuron, exclneurons) make NaN and fix
%                     code
                    frse(ineuron) = sum(firings(:,2) == ineuron) / (cfg.Nt * cfg.dtsc / 1000);
                end
                [~, basedensity, frexcitatory] = MarijnEasyHist(cfg,firings,icell,exclneurons);
                ii        = 1:length(basedensity) - 1;
                iip       = find(basedensity(ii) <  cfg.PeakThr & basedensity(ii + 1) >= cfg.PeakThr);
                iim       = find(basedensity(ii) >= cfg.PeakThr & basedensity(ii + 1) <  cfg.PeakThr);
                len       = min([length(iim) length(iip)]);
                if len > 0 
                    iip       = iip(1:len,1);
                    iim       = iim(1:len,1);
                end
                numburst  = length(iim);
                range     = nan(numburst,2);
                baseamps  = nan(numburst,1);
                burstrate = numburst / (length(basedensity) * cfg.dtsc / 1000);
                for i = 1:numburst
                    range(i,:)  = [iip(i) iim(i)];
                    baseamps(i) = max(basedensity(range(i,1):range(i,2)));
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%
            % stimulation responses %
            %%%%%%%%%%%%%%%%%%%%%%%%%
            stimamps = NaN; 
            stimlats = NaN; 
            TP       = zeros(length(cfg.PeakThrROC),1); 
            FP       = zeros(length(cfg.PeakThrROC),1);
            if cfg.patdrive > 0
                stimfirings = RStore(inettyp, inet, ipara).stimfirings;
                [~, stimdensity,~] = MarijnEasyHist(cfg,stimfirings,icell,exclneurons);
                stimon = find(diff(cfg.Ipatdrive) == 1);
                stimon = stimon(stimon + cfg.stimwindow(2) / cfg.dtsc < cfg.Nt);
                for istim = 1:length(stimon)
                    stimtims = stimon(istim) + cfg.stimwindow(1) / cfg.dtsc:stimon(istim) + cfg.stimwindow(2) / cfg.dtsc;
                    if cfg.needsbasefire
                        [stimbaseamps(istim) ~]    	  = max(basedensity(stimtims));
                    end
                    [stimamps(istim) stimlats(istim)] = max(stimdensity(stimtims));
                    if stimamps(istim) < cfg.PeakThr && cfg.patdrive == 1 % only calc latency/amplitude for suc bursts
                        stimamps(istim) = NaN;
                        stimlats(istim) = NaN;
                    end
                end
                stimlats = stimlats .* cfg.dtsc;
                if cfg.needsbasefire
                    for ithr = 1:length(cfg.PeakThrROC)
                        for istim = 1:length(stimon)
                            if stimbaseamps(istim) > cfg.PeakThrROC(ithr)
                                FP(ithr) = FP(ithr) + 1 / (length(stimon));
                            end
                            if stimamps(istim) > cfg.PeakThrROC(ithr)
                                TP(ithr) = TP(ithr) + 1 / (length(stimon));
                            end
                        end
                    end
                end
            end
            SStore(inettyp,inet,ipara).fr           = frexcitatory;
            SStore(inettyp,inet,ipara).frse         = frse;
            SStore(inettyp,inet,ipara).burstrate    = burstrate;
            SStore(inettyp,inet,ipara).range        = range;
            SStore(inettyp,inet,ipara).baseamps     = baseamps;
            SStore(inettyp,inet,ipara).stimamps     = stimamps;
            SStore(inettyp,inet,ipara).stimamp      = stimamps(1);
            SStore(inettyp,inet,ipara).stimlat      = stimlats(1);
            SStore(inettyp,inet,ipara).tp           = TP;
            SStore(inettyp,inet,ipara).fp           = FP;
        end
    end
end %#ok<*AGROW>