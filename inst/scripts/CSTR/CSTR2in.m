function [Fvec, CA0vec, T0vec, Tcinvec, Fcvec] = ...
                 CSTR2in(t, condition, tau)

%  Last modified 2 May 2005

n       = length(t);
CA0vec  =   2*ones(n,1);
T0vec   = 323*ones(n,1);
Fcvec   =  15*ones(n,1);

if nargin < 3, tau = 1;  end

switch condition
    
    case 'all_cool_step'
        
        %  compute F
        
        Fvec = 1.0*ones(n,1);
        Fvec(find(4 <= t & t <  8)) = 1.5;
        Fvec(find(8 <= t & t < 12)) = 0.5;
        
        %  compute C_{A0}
        
        CA0vec = 2.0*ones(n,1);
        CA0vec(find(16 <= t & t < 20)) = 2.2;
        CA0vec(find(20 <= t & t < 24)) = 1.8;
        
        %  compute T0
        
        T0vec = 323*ones(n,1);
        T0vec(find(28 <= t & t < 32)) = 343;
        T0vec(find(32 <= t & t < 36)) = 303;
        
        %  compute Tcin
        
        Tcinbase = 335;
        Tcinvec  = Tcinbase*ones(n,1);
        Tcinvec(find(40 <= t & t < 44)) = Tcinbase+5;
        Tcinvec(find(44 <= t & t < 48)) = Tcinbase-5;
        
        %  compute Fc
        
        Fcvec = 15*ones(n,1);
        Fcvec(find(52 <= t & t < 56)) = 20;
        Fcvec(find(56 <= t & t < 60)) = 10;
        
    case 'all_hot_step'
        
        %  compute F
        
        Fvec = 1.0*ones(n,1);
        Fvec(find(4 <= t & t <  8)) = 1.5;
        Fvec(find(8 <= t & t < 12)) = 0.5;
        
        %  compute C_{A0}
        
        CA0vec = 2.0*ones(n,1);
        CA0vec(find(16 <= t & t < 20)) = 2.2;
        CA0vec(find(20 <= t & t < 24)) = 1.8;
        
        %  compute T0
        
        T0vec = 323*ones(n,1);
        T0vec(find(28 <= t & t < 32)) = 343;
        T0vec(find(32 <= t & t < 36)) = 303;
        
        %  compute Tcin
        
        Tcinbase = 365;
        Tcinvec  = Tcinbase*ones(n,1);
        Tcinvec(find(40 <= t & t < 44)) = Tcinbase+5;
        Tcinvec(find(44 <= t & t < 48)) = Tcinbase-5;
        
        %  compute Fc
        
        Fcvec = 15*ones(n,1);
        Fcvec(find(52 <= t & t < 56)) = 20;
        Fcvec(find(56 <= t & t < 60)) = 10;
        
    case 'all_hot_ramp'

        Fvec    = 1.0*ones(n,1);
        Tcinvec = 365*ones(n,1);
        
        index = find( 2 <= t & t < 10);
        Fvec(index) = 1.5;
        index = find(t >= 10 & t < 14);
        Fvec(index) = -0.2*(t(index)-10) + 1.5;
        index = find(t >= 14 & t < 18);
        Fvec(index) = -0.2*4 + 1.5;
        
        index = find(t >= 26 & t < 34);
        CA0vec(index) = 2.5;
        index = find(t >= 34 & t < 38);
        CA0vec(index) = -0.2*(t(index)-34) + 2.5;
        index = find(t >= 38 & t < 42);
        CA0vec(index) = -0.2*4 + 2.5;

        index = find(t >= 50 & t < 58);
        T0vec(index) = 353;
        index = find(t >= 58 & t < 62);
        T0vec(index) = -20*(t(index)-58) + 353;
        index = find(t >= 62 & t < 66);
        T0vec(index) = -20*4 + 353;
        
        index = find(t >= 74 & t < 82);
        Tctop = 390;
        Tcinvec(index) = Tctop;
        index = find(t >= 82 & t < 86);
        Tcinvec(index) = -10*(t(index)-82) + Tctop;
        index = find(t >= 86 & t < 90);
        Tcinvec(index) = -10*4 + Tctop;
        
        index = find(t >=  98 & t < 106);
        Fcvec(index) = 25;
        index = find(t >= 106 & t < 110);
        Fcvec(index) = -5*(t(index)-106) + 25;
        index = find(t >= 110 & t < 114);
        Fcvec(index) = -5*4 + 25;
        
    case 'all_cool_ramp'

        Fvec    = 0.05*ones(n,1);
        Tcinvec = 330*ones(n,1);
        
        index = find( 2 <= t & t < 10);
        Fvec(index) = 1.5;
        index = find(t >= 10 & t < 14);
        Fvec(index) = -0.2*(t(index)-10) + 1.5;
        index = find(t >= 14 & t < 18);
        Fvec(index) = -0.2*4 + 1.5;
        
        index = find(t >= 26 & t < 34);
        CA0vec(index) = 2.5;
        index = find(t >= 34 & t < 38);
        CA0vec(index) = -0.2*(t(index)-34) + 2.5;
        index = find(t >= 38 & t < 42);
        CA0vec(index) = -0.2*4 + 2.5;

        index = find(t >= 50 & t < 58);
        T0vec(index) = 353;
        index = find(t >= 58 & t < 62);
        T0vec(index) = -20*(t(index)-58) + 353;
        index = find(t >= 62 & t < 66);
        T0vec(index) = -20*4 + 353;
        
        index = find(t >= 74 & t < 82);
        Tctop = 355;
        Tcinvec(index) = Tctop;
        index = find(t >= 82 & t < 86);
        Tcinvec(index) = -10*(t(index)-82) + Tctop;
        index = find(t >= 86 & t < 90);
        Tcinvec(index) = -10*4 + Tctop;
        
        index = find(t >=  98 & t < 106);
        Fcvec(index) = 25;
        index = find(t >= 106 & t < 110);
        Fcvec(index) = -5*(t(index)-106) + 25;
        index = find(t >= 110 & t < 114);
        Fcvec(index) = -5*4 + 25;
        
    case 'Tc_hot_exponential'
        
        Fvec    = 1.0*ones(n,1);
        Tcinvec = 365*ones(n,1);
        
        index = find(t < 10); 
        Tcinvec(index) = 400 - (400 - 365).*exp(-(t(index)     )./tau);
        index = find(10 <= t & t < 20); 
        Tcinvec(index) = 344 - (344 - 400).*exp(-(t(index) - 10)./tau);
        index = find(20 <= t); 
        Tcinvec(index) = 365 - (365 - 344).*exp(-(t(index) - 20)./tau);

    case 'Tc_cool_exponential'
        
        Fvec    = 0.05*ones(n,1);
        Tcinvec =  330*ones(n,1);

        index = find(t/5 < 10); 
        Tcinvec(index) = 400 - (400 - 365).*exp(-(t(index)     )./tau);
        index = find(10 <= t/5 & t/5 < 20); 
        Tcinvec(index) = 344 - (344 - 400).*exp(-(t(index) - 10)./tau);
        index = find(20 <= t/5); 
        Tcinvec(index) = 365 - (365 - 344).*exp(-(t(index) - 20)./tau);

     case 'Tc_hot_ramp'

        Fvec    = 1.0*ones(n,1);
        Tcinvec = 365*ones(n,1);
        
        index = find(2 <= t & t < 10); 
        Tcinvec(index) = 400;
        index = find(10 <= t & t < 14); 
        Tcinvec(index) = -14*(t(index)-10) + 400;
        index = find(14 <= t & t < 18); 
        Tcinvec(index) = -14*4 + 400;
        
     case 'Tc_cool_ramp'

        Fvec    = 0.05*ones(n,1);
        Tcinvec =  330*ones(n,1);
        
        index = find(2 <= t/5 & t/5 < 10); 
        Tcinvec(index) = 340;
        index = find(10 <= t/5 & t/5 < 14); 
        Tcinvec(index) = -0.8*(t(index)-5*10) + 340;
        index = find(14 <= t/5 & t/5 < 22); 
        Tcinvec(index) = -5*0.8*4 + 340;
        
     case 'Tc_hot_step'

        Fvec    = 1.0*ones(n,1);
        Tcinvec = 365*ones(n,1);
        
        index = find(2 <= t & t < 12); 
        Tcinvec(index) = 350;
        
     case 'Tc_cool_step'

        Fvec    = 0.05*ones(n,1);
        Tcinvec =  335*ones(n,1);
        
        index = find(2 <= t & t < 12); 
        Tcinvec(index) = 320;
        
end

