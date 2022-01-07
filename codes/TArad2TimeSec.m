function Result=TArad2TimeSec(AngMom, Ecc, TArad, mu)

switch 1
    case Ecc < 1  % Elliptic Orbit, a > 0
        aSemMaj=AngMom^2/mu/(1-Ecc^2);
        Period=2*pi*sqrt(aSemMaj^3/mu);
        EccAnomaly=2*atan( sqrt((1-Ecc)/(1+Ecc)) * tan(TArad/2) ); 
        MeanAnomaly=EccAnomaly-Ecc*sin(EccAnomaly);
        Result=Period/(2*pi)*MeanAnomaly;
        if Result < 0
            Result=Period+Result;
        end;
    case Ecc == 1 % Parabolic Trajectory, a -> Inf
        TanTA2=tan(TArad/2);
        Result=AngMom^3/mu^2*( TanTA2/2 + TanTA2^3/6 );
    case Ecc > 1  % Hyperbolic Trajectory, a < 0
        HypEccAnomaly  = 2*atanh( sqrt((Ecc-1)/(Ecc+1)) * tan(TArad/2) );
        HypMeanAnomaly = Ecc*sinh(HypEccAnomaly)-HypEccAnomaly;
        Result=AngMom^3/mu^2*(Ecc^2-1)^(-1.5)*HypMeanAnomaly;   
end;

