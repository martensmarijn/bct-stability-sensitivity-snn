function cfg = MarijnDefaultCfg
% default network structure
cfg.Ne           = 480;
cfg.normweights  = false;
cfg.iparidx      = [1 1]; % set which index in matrix changed for cfg.var
cfg.SynConMat    = [.013   0.013  0.007;
                    .010   0.013  0.004;
                    .004   0      0];
cfg.pconSyn      = [.05     1      1; % (0.05 PC-PC) http://jp.physoc.org/content/551/1/139.abstract, (1 PV PC): http://www.jneurosci.org/content/31/37/13260.full and (1 SST PC): http://www.jneurosci.org/content/31/37/13260.full#ref-16 and scanziani paper
                    .5      1      1; % (1 PC-PV and PC-SST are reciprocal) http://jp.physoc.org/content/551/1/139.abstract
                    .5      0      0];% IMPORTANT: inhibitory connectivity is always fully connected, see MarijnGenNetwork to implement partial connectivity.

% default neuron params
cfg.tAMPA       = 2;
cfg.tGABA       = 10;
cfg.EGABA       = -80;
cfg.EAMPA       = 0;

% default plasticity params
cfg.taustdp    = 20;
cfg.homtau      = 2;
cfg.Pamp        = 0.005;
cfg.Damp        = 0.00525;
cfg.dostdp      = false;
cfg.stdpselect_t= 1E9;
cfg.frtar       = 1;
cfg.morrisonstdp = false;
cfg.w0fac       = 80;

% default stimulation params
cfg.patdrive    = 0;
cfg.NFLEX_stim  = 6;
cfg.CurFlex     = 8;
cfg.StimDuration= 25;
cfg.IstimAmpRandom = 15;
cfg.Ibg         = [3 0 0];
cfg.Isd         = [.5 .5 .5]; %std in noise current
cfg.stimwindow  = [5 65];
cfg.noisehz     = 0.07; % random spikes per neuron (Hz)
cfg.needsbasefire = true;
cfg.only1basefire = false;

% default simulation params
cfg.loadwhattype = 0;
cfg.loadwhatipara = 0;
cfg.dtsc        = 0.1;
cfg.Nt          = 2500;
cfg.npara       = 10;
cfg.inettypes   = 1:4;
cfg.Nnet        = 1;
cfg.vary_spike_with_inet = 3;
cfg.PeakThr     = 10;
cfg.PeakThrROC  = [linspace(1,20,39) 25 30 40 50 100 150 200];

% default plot params
cfg.plot.times        = 5000:12000;
cfg.plot.sig_in_ms2   = 25;
cfg.plot.bint         = .5;
cfg.plot.hp.LineWidth = 2;
cfg.plot.hp.FontSize  = 8;
cfg.plot.h.linewidth  = 2;
cfg.plot.titlar       = {'PCOR','ACOR','UCOR','XCOR1','XCOR2','XCOR3'};
cfg.plot.clr          = {'b','g','r','c','k','m'};
cfg.plot.clrrgb       = [0 0 1; 0 1 0; 1 0 0; 0 1 1; 0 0 0; 1 0 1];
cfg.plot.nrnclr       = {'k','r','c'};

% default netstats extract params
cfg.do.long      = false;
cfg.do.path      = false;
cfg.do.kcore     = false;
cfg.do.distance  = false;
cfg.do.kcorecentr= false;
cfg.qmax         = 3;
cfg.kdegmin      = 34;
end
