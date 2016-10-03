addpath(genpath('/home/marijn/documents/git/bct'))
datdir = '/home/marijn/documents/bct/data';
figdir ='/home/marijn/documents/bct/figures';

for whattype = [1]
whattype
tic
% try
    clearvars -except *dir whattype
    GenNets           = 1;
    ExtractNetStats   = 1;
    ExtractSpikeStats = 1;
    StoreData         = 1;

    cfg          = MarijnDefaultCfg;
    cfg.whattype = whattype;
    cfg          = MarijnCaseCfg(cfg);

    if GenNets
        gcp; %create matlabpool
        RStoreGen = cell(cfg.Nnet,1);
        parfor inet = 1:cfg.Nnet
%         for inet = 1:cfg.Nnet
            inet
            cfg          = MarijnDefaultCfg;
            cfg.whattype = whattype;
            cfg          = MarijnCaseCfg(cfg);
            MStore       = [];
            RStore       = [];
            for inettyp = cfg.inettypes
                con          = MarijnLoadCon(cfg,inettyp,inet,datdir);
                for ipara = 1:cfg.npara
                    cfg    = MarijnBuildFinalCfg(cfg, ipara);
                    MStore = MarijnGenNetwork(cfg, inettyp, inet, ipara, MStore);
                    cfg    = MarijnGenStimulus(cfg);
                    RStore = MarijnSpiking(cfg, inettyp, inet, ipara, MStore, RStore, datdir);
                end
            end
            RStoreGen{inet} = RStore;
        end
        clearvars -except *dir whattype RStoreGen *Data Extract*
        RStore       = [];
        MStore       = [];
        cfg          = MarijnDefaultCfg;
        cfg.whattype = whattype;
        cfg          = MarijnCaseCfg(cfg);
        for inet = 1:cfg.Nnet
            for inettyp = cfg.inettypes
                for ipara = 1:cfg.npara      
                    cfg    = MarijnBuildFinalCfg(cfg, ipara);
                    MStore = MarijnGenNetwork(cfg, inettyp, inet, ipara, MStore);
                    cfg    = MarijnGenStimulus(cfg);
                    if cfg.needsbasefire
                        RStore(inettyp,inet,ipara).firings = RStoreGen{inet}(inettyp,ipara).firings;
                    end
                    if cfg.patdrive > 0
                        RStore(inettyp,inet,ipara).stimfirings = RStoreGen{inet}(inettyp,ipara).stimfirings;
                    end
                    if cfg.dostdp
                        RStore(inettyp,inet,ipara).con   = RStoreGen{inet}(inettyp,ipara).con;
                        RStore(inettyp,inet,ipara).constdp = RStoreGen{inet}(inettyp,ipara).constdp;
                    end
                end
            end
            clearvars -except *dir whattype RStoreGen *Data Extract* cfg RStore MStore
        end
    else
        load(sprintf('%s/MarijnSpiking_%d.mat',datdir,cfg.whattype));
    end
%% storing data and extracting stats
    if StoreData
        if ~exist('SStore','var'); SStore = []; end
        eval(sprintf('save %s/MarijnSpiking_%d cfg MStore RStore SStore -v7.3',datdir,cfg.whattype));
    end
%     cfg          = MarijnDefaultCfg;
%     cfg.whattype = whattype;
%     cfg          = MarijnCaseCfg(cfg); 
    if ExtractNetStats
        MStore = MarijnExtractNetStats(cfg, MStore);
    end
    if ExtractSpikeStats
        SStore = MarijnExtractSpikeStats(cfg, RStore, MStore);
    end
    if StoreData
        eval(sprintf('save %s/MarijnSpiking_%d cfg MStore RStore SStore -v7.3',datdir,cfg.whattype));
    end
%% plotting
    for optidx = 1:length(cfg.plot.options)
        cfg.plot.icf(optidx) = figure;
        for inettyp = cfg.inettypes
            for inet = 1:cfg.Nnet
                for ipara = 1:cfg.npara
                    cfg = MarijnBuildFinalCfg(cfg, ipara);
                    cfg = MarijnGenStimulus(cfg);
                    MarijnPlotParaStuff(cfg, inettyp, inet, ipara, MStore, RStore, SStore, optidx);
                end
                MarijnPlotNetStuff(cfg, inettyp, inet, MStore, RStore, SStore, optidx);
            end
        end
        MarijnPlotMetaStuff(cfg, MStore, SStore, RStore, optidx);
        eval(sprintf('print -dpdf %s/MarijnSpiking_%d_%s',figdir,cfg.whattype,cfg.plot.options{optidx}));
%         export_fig('filename',[sprintf('%s/MarijnSpiking_%d_%s',figdir,cfg.whattype,cfg.plot.options{optidx}) '.pdf'],'-pdf')
    end
toc
% catch
%     toc
%     lasterror
% end
end
