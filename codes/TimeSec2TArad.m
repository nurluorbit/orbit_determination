function Result=TimeSec2TArad(AngMom, Ecc, TimeSec, mu, Precision)

switch 1
    case Ecc < 1  % Elliptic Orbit, a > 0
        aSemMaj=AngMom^2/mu/(1-Ecc^2);
        Period=2*pi*sqrt(aSemMaj^3/mu);
        if TimeSec > Period;
            TimeSec = TimeSec-Period*fix(TimeSec/Period);
        end
        MeanAnomaly=2*pi/Period*TimeSec;
           EccAnomalyOld=MeanAnomaly;
           EccAnomalyNew=MeanAnomaly+100;
           while abs(EccAnomalyOld-EccAnomalyNew)>Precision
               EccAnomalyOld=EccAnomalyNew;
               EccAnomalyNew=MeanAnomaly+Ecc*sin(EccAnomalyOld);
           end;
        Result1=2*atan( sqrt( (1+Ecc)/(1-Ecc) ) *tan(EccAnomalyNew/2) );
        Result=Result1;
        
        if Result1 < 0
            Result=Result1+2*pi;
        end;
        
    case Ecc == 1 % Parabolic Trajectory, a -> Inf
        ParMeanAnomaly=mu^2/AngMom^3*TimeSec;
        ThreePMA=3*ParMeanAnomaly;
        Result1=(ThreePMA+sqrt(ThreePMA^2+1))^(1/3);
        Result=2*atan(Result1-1/Result1);
        
    case Ecc > 1  % Hyperbolic Trajectory, a < 0
        HypMeanAnomaly=mu^2/AngMom^3*(Ecc^2-1)^(1.5)*TimeSec;
        
           HypEccAnomalyOld=HypMeanAnomaly;
           HypEccAnomalyNew=HypMeanAnomaly+100;
           while abs(HypEccAnomalyOld-HypEccAnomalyNew)>Precision
               HypEccAnomalyOld=HypEccAnomalyNew;
               HypEccAnomalyNew=asinh((HypEccAnomalyOld+HypMeanAnomaly)/Ecc);
           end;
        Result=2*atan( sqrt( (Ecc+1)/(Ecc-1) )*tanh(HypEccAnomalyNew/2) );
end;

