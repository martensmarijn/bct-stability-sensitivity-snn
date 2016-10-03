function cfg = MarijnGenStimulus(cfg)

Ipatdrive   = zeros(cfg.Nt, 1);
            
if cfg.patdrive > 0
    switch cfg.patdrive
        case{1}
            iit             = 150:1000:5E4;
            iid             = round(cfg.StimDuration / cfg.dtsc);
            ii              = find(iit / cfg.dtsc + iid < cfg.Nt);
            for i=1:length(ii)
                range            = round(iit(ii(i)) / cfg.dtsc) + (1:iid);
                Ipatdrive(range) = 1;
            end
            Ipatdrive       = Ipatdrive(1:cfg.Nt);
        case{2}
            iit             = 300:70:5E4;
            iid             = round(cfg.StimDuration / cfg.dtsc);
            ii              = find(iit / cfg.dtsc + iid < cfg.Nt);
            for i=1:length(ii)
                range            = round(iit(ii(i)) / cfg.dtsc) + (1:iid);
                Ipatdrive(range) = 1;
            end
            Ipatdrive       = Ipatdrive(1:cfg.Nt);    
    	case{3}
            iit             = 300:200:5E4;
            iid             = round(cfg.StimDuration / cfg.dtsc);
            ii              = find(iit / cfg.dtsc + iid < cfg.Nt);
            for i=1:length(ii)
                range            = round(iit(ii(i)) / cfg.dtsc) + (1:iid);
                Ipatdrive(range) = 1;
            end
            Ipatdrive       = Ipatdrive(1:cfg.Nt); 
         case{4}
            iit             = 300:500:10E4;
            iid             = round(cfg.StimDuration / cfg.dtsc);
            ii              = find(iit / cfg.dtsc + iid < cfg.Nt);
            for i=1:length(ii)
                range            = round(iit(ii(i)) / cfg.dtsc) + (1:iid);
                Ipatdrive(range) = 1;
            end
            Ipatdrive       = Ipatdrive(1:cfg.Nt);    
            
    end
end

cfg.Ipatdrive   = Ipatdrive;
