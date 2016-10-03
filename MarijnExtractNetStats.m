function MStore = MarijnExtractNetStats(cfg, MStore)    

for inettyp = cfg.inettypes
    for inet = 1:size(MStore,2)
        for ipara = 1:size(MStore, 3)
            cfg = MarijnBuildFinalCfg(cfg, ipara);
            M   = full(MStore(inettyp,inet,cfg.netpara).M(1:cfg.Ne,1:cfg.Ne));
            if cfg.do.distance
                disp('doing distance calcs...')
                SelectOffDiag = ones(cfg.Ne);
                SelectOffDiag = SelectOffDiag - diag(diag(SelectOffDiag));
                ii            = find(SelectOffDiag > 0); %linear indices;
                D             = distance_bin(M);
                Dlin          = D(ii);
                Dlin          = Dlin(find(isfinite(Dlin)));
                Dcount        = full(sparse(Dlin,1,1));
                MeanD         = mean(Dlin);
                MStore(inettyp,inet,ipara).Dcount    	= Dcount;
                MStore(inettyp,inet,ipara).MeanD      	= MeanD;
            end
            if cfg.do.kcore
                disp('doing kcore calcs...')
                [CIJkcore, kn, peelorder, peellevel] 	= kcore_bd(M, cfg.kdegmin);
                MStore(inettyp,inet,ipara).kn        	= kn;
            end
            if cfg.do.kcorecentr
                disp('doing kcore centrality calcs...')
                [coreness,kncentr]                      = kcoreness_centrality_bu(M);
                MStore(inettyp,inet,ipara).kncentr     	= kncentr;
            end
            if cfg.do.long 
                disp('doing edge betweenness calcs...')
                [EBC BC] = edge_betweenness_bin(M);
                MStore(inettyp,inet,ipara).BC    = BC;
                MStore(inettyp,inet,ipara).MaxBC = max(BC);
                MStore(inettyp,inet,ipara).MinBC = min(BC);
                % [EC,ec,degij] = edge_nei_overlap_bd(M); %I wonder what you can do with this
            end
            if cfg.do.path
                InD                 = MStore(inettyp,inet,cfg.netpara).InDegree;
                OD                  = MStore(inettyp,inet,cfg.netpara).OutDegree;
                [ODsorted IndxOD]   = sort(OD,'Ascend'); 
                [InDsorted IndxInD] = sort(InD,'Ascend');
%                 sources           = [IndxOD(1:cfg.HowManyStim)]; % setting was cfg.HowManyStim  = 10;
                sources             = [IndxOD(1) IndxOD(round(cfg.Ne/2)) IndxOD(cfg.Ne) IndxInD(1) IndxInD(round(cfg.Ne/2)) IndxInD(cfg.Ne)]  ;
                [Pq,tpath,plq,qstop,allpths,util] = findpaths(M,sources,cfg.qmax,0);
                Pq                  = Pq(sources,:,:);
                disp('doing path finding calcs...')
                MStore(inettyp,inet,ipara).sources = sources;
                %MStore(inettyp,inet,ipara).Pq      = Pq; %perhaps better not  to save
                MStore(inettyp,inet,ipara).Pqs1    = squeeze(sum(Pq>0,2));
                MStore(inettyp,inet,ipara).Pqs2    = squeeze(sum(Pq,2));
            end
        end
    end
end
