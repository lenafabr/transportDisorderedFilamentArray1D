function welldepth = getWellDepthPrediction(params,MTparam,brownbridge,lowestorder)
% compute the well depth prediction from a set of kinetic and
% configurational parameters

beta = params.v/params.kd;
gamma = sqrt(params.D/(MTparam.ncross*params.ka));
Ldom = MTparam.domain(end)-MTparam.domain(1);
Lmt = MTparam.avglen;
ncross = MTparam.ncross;

if (brownbridge)
    welldepth = sqrt(pi/2)*sqrt(2*Ldom*Lmt/ncross)*beta/(beta^2+gamma^2);
else
    welldepth = sqrt(8/pi)*sqrt(2*Ldom*Lmt/ncross)*beta/(beta^2+gamma^2);
end

if (lowestorder)
    nterm = sqrt(1/3);
else
    nterm = sqrt(1/3 + 1/2/ncross + 1/4/ncross^2)
end

welldepth = welldepth*nterm;

end