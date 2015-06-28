function [fM, WM, TM] = getConstPBattMap(Pbatt, WeM, Wreq, Treq, enginemap_trq,enginemap_spd,enginemap)

    nm = 0.85;
    ng = 0.85;
    S = 30;
    Rr = 78;
    R = Rr/S;
    K = 4.113;
    kcom = [1 1 -1 -1;1 -1 1 -1];
    
    fM = nan(numel(WeM), 1);
    WM = nan(numel(WeM), 1);
    TM = nan(numel(WeM), 1);
    for k = 1:length(WeM)

        We = WeM(k);
        Wm = K*Wreq;
        Wg = -R*K*Wreq + (1+R)*We;

        for l = 1:length(kcom)

            km = kcom(1,l);
            kg = kcom(2,l);
            mtx = [ng^kg*Wg nm^km*Wm 0; 0 1 R/(1+R); 1 0 1/(1+R)];
            T = inv(mtx)*[Pbatt; Treq/K; 0];

            % verify this!
            Tgo = T(1);
            Tmo = T(2);
            Teo = T(3);

            % find T*W that correlates the the right value of k (-1 or 1)

            % Problem 1 : more than 1 combination of km & kg that work
            % Answer  1 : Check cover

            % Problem 2 : If no Tm & Tg fullfil any of the condition, it will remain the same as the Previous Tm & Tg
            % Answer  2 : Check cno

            if (Tmo*Wm*km <= 0)&&(Tgo*Wg*kg <= 0)
                Tm = Tmo;
                Tg = Tgo;
                Te = Teo;
%                 countin = countin + 1;
%                         fprintf('%d %d %d, km = %d, kg = %d\n', i, j, k, km, kg)
                if Tg > 0;
                    debug = 1;
                end
            else
%                 countout = countout + 1;
            end

            % Check if NO Combination works!
%             if countout == 4
%                 c = c + 1;
%                 cno(c,:) = [i j k];
                % Check if more than 1 Combinations work!
%             elseif countin > 1
%                 d = d + 1;
%                 cover(d,:) = [i j k];
%             end

        end

%         countin = 0;
%         countout = 0;
        if Wreq > 0
            debug = 1;
        end

        fM(k) = interp2(enginemap_trq,enginemap_spd,enginemap,Te,We);
        TM(k) = Te;
        WM(k) = We;
    end