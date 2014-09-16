%CalibImpact
%The Climate Framework for Uncertainty, Negotiation and Distribution,
%version 4.0-matlab-national
%
%This script is part of FUND 4.0 NG
%It calibrates five alternative aggregrate damage functions
%
%Richard Tol, 9 September 2014
%This code is protected by the MIT License

%global impact
vtemp = observedImpact(:,1);
vimp = observedImpact(:,2);
vN = length(vtemp);
vtemp = vtemp(1:vN-1);
vimp = vimp(1:vN-1);

imppar = zeros(NImpact,3);

[imppar(1,1) imppar(1,2) imppar(1,3)] = fittol(vtemp,vimp);
[imppar(2,1) imppar(2,2) imppar(2,3)] = fitweitzman(vtemp,vimp);
[imppar(3,1) imppar(3,3)] = fitnordhaus(vtemp,vimp);
[imppar(4,1) imppar(4,3)] = fithope(vtemp,vimp);
[imppar(5,1) imppar(5,3)] = fitploeg(vtemp,vimp);

%national impact
ctrimppar = zeros(NImpact,NCountry,2);
xnatimp = zeros(NImpact,NCountry);

for c = 1:NCountry,
    ctrimppar(1,c,1) = imppar(1,1)*(YpCCtr(c,NHistYear,1)/YpC2010)^impincelas;
    ctrimppar(1,c,2) = imppar(1,2)*(YpCCtr(c,NHistYear,1)/YpC2010)^impincelas;
    xnatimp(1,c) = ctrimppar(1,c,1)*2.5 + ctrimppar(1,c,2)*2.5^2;
    ctrimppar(2,c,1) = imppar(2,1)*(YpCCtr(c,NHistYear,1)/YpC2010)^impincelas;
    ctrimppar(2,c,2) = imppar(2,2)*(YpCCtr(c,NHistYear,1)/YpC2010)^impincelas;
    xnatimp(2,c) = ctrimppar(2,c,1)*2.5^2 + ctrimppar(2,c,2)*2.5^6;
    ctrimppar(3,c,1) = imppar(3,1)*(YpCCtr(c,NHistYear,1)/YpC2010)^impincelas;
    xnatimp(3,c) = ctrimppar(3,c,1)*2.5^2;
    ctrimppar(4,c,1) = imppar(4,1)*(YpCCtr(c,NHistYear,1)/YpC2010)^impincelas;
    xnatimp(4,c) = ctrimppar(4,c,1)*2.5^1.3;
    ctrimppar(5,c,1) = imppar(5,1)*(YpCCtr(c,NHistYear,1)/YpC2010)^impincelas;
    xnatimp(5,c) = ctrimppar(5,c,1)*(exp(2*2.5/4.33)-1);
end

%rescale to match global impact for 2010 income and 2.5K warming
xtotimp = sum(xnatimp(1,:)'.*YCtr(:,NHistYear,1));
xglbimp = (imppar(1,1)*2.5 + imppar(1,2)*(2.5^2))*Y(NHistYear,1);
ctrimppar(1,:,1) = ctrimppar(1,:,1)*xtotimp/xglbimp;
ctrimppar(1,:,2) = ctrimppar(1,:,2)*xtotimp/xglbimp;

xtotimp = sum(xnatimp(2,:)'.*YCtr(:,NHistYear,1));
xglbimp = (imppar(2,1)*(2.5^2) + imppar(2,2)*(2.5^6))*Y(NHistYear,1);
ctrimppar(2,:,1) = ctrimppar(2,:,1)*xtotimp/xglbimp;
ctrimppar(2,:,2) = ctrimppar(2,:,2)*xtotimp/xglbimp;

xtotimp = sum(xnatimp(3,:)'.*YCtr(:,NHistYear,1));
xglbimp =imppar(3,1)*(2.5^2)*Y(NHistYear,1);
ctrimppar(3,:,1) = ctrimppar(3,:,1)*xtotimp/xglbimp;

xtotimp = sum(xnatimp(4,:)'.*YCtr(:,NHistYear,1));
xglbimp =imppar(4,1)*(2.5^1.3)*Y(NHistYear,1);
ctrimppar(4,:,1) = ctrimppar(4,:,1)*xtotimp/xglbimp;

xtotimp = sum(xnatimp(5,:)'.*YCtr(:,NHistYear,1));
xglbimp =imppar(5,1)*(exp(2*2.5/4.33)-1)*Y(NHistYear,1);
ctrimppar(5,:,1) = ctrimppar(5,:,1)*xtotimp/xglbimp;

clear x*