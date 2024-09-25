function [ wDenom ] = addWienerDenominator_3D(wd,sp,d,b)
%set up the Wiener Filter 
%   此处显示详细说明
siz=size(wd);
w=siz(2);
h=siz(1);
dir=sp.Dir(d);
cyclMicron=sp.cyclesPerMicron;
cnt=[h/2+1,w/2+1];

x=1:w;
y=1:h;
[x,y]=meshgrid(x,y);

rad1=hypot(x-cnt(2)-(b-1)*dir.px,y-cnt(1)-(b-1)*dir.py)*cyclMicron;
rad2=hypot(x-cnt(2)+(b-1)*dir.px,y-cnt(1)+(b-1)*dir.py)*cyclMicron;


otfVal1=abs(getOtfVal(sp.OtfProvider,rad1)).^2;
otfVal2=abs(getOtfVal(sp.OtfProvider,rad2)).^2;


if b==1 
%     otfVal1=otfVal1.*valAttenuation(rad1,sp.OtfProvider.attStrength/1.05,1.0*sp.OtfProvider.attFWHM)/2;
%     otfVal2=otfVal2.*valAttenuation(rad1,sp.OtfProvider.attStrength/1.05,1.0*sp.OtfProvider.attFWHM)/2;
    otfVal1=otfVal1.*valAttenuation(rad1,sp.OtfProvider.attStrength/1.0,1.0*sp.OtfProvider.attFWHM)/2;
    otfVal2=otfVal2.*valAttenuation(rad1,sp.OtfProvider.attStrength/1.0,1.0*sp.OtfProvider.attFWHM)/2;
elseif b==2
%     otfVal1=1*otfVal1.*valAttenuation(rad1,sp.OtfProvider.attStrength/1.0,1*sp.OtfProvider.attFWHM)/1.0;
%     otfVal2=1*otfVal2.*valAttenuation(rad1,sp.OtfProvider.attStrength/1.0,1*sp.OtfProvider.attFWHM)/1.0;
%     otfVal1=otfVal1.*valAttenuation(rad1,sp.OtfProvider.attStrength/1.15,0.85*sp.OtfProvider.attFWHM)/1.2;
%     otfVal2=otfVal2.*valAttenuation(rad1,sp.OtfProvider.attStrength/1.15,0.85*sp.OtfProvider.attFWHM)/1.2;

    otfVal1=otfVal1.*valAttenuation(rad1,sp.OtfProvider.attStrength/1.0,1.0*sp.OtfProvider.attFWHM)/1.2;
    otfVal2=otfVal2.*valAttenuation(rad1,sp.OtfProvider.attStrength/1.0,1.0*sp.OtfProvider.attFWHM)/1.2;
else
%     otfVal1=1*otfVal1.*valAttenuation(rad1,sp.OtfProvider.attStrength/1.0,1*sp.OtfProvider.attFWHM)/1.0;
%     otfVal2=1*otfVal2.*valAttenuation(rad1,sp.OtfProvider.attStrength/1.0,1*sp.OtfProvider.attFWHM)/1.0;

%   otfVal1=1*otfVal1.*valAttenuation(rad1,sp.OtfProvider.attStrength/1.15,0.85*sp.OtfProvider.attFWHM)/1.2;
%   otfVal2=1*otfVal2.*valAttenuation(rad1,sp.OtfProvider.attStrength/1.15,0.85*sp.OtfProvider.attFWHM)/1.2; 
    
    otfVal1=1*otfVal1.*valAttenuation(rad1,sp.OtfProvider.attStrength/1.0,1.0*sp.OtfProvider.attFWHM)/1.2;
    otfVal2=1*otfVal2.*valAttenuation(rad1,sp.OtfProvider.attStrength/1.0,1.0*sp.OtfProvider.attFWHM)/1.2; 
end

wd=wd+otfVal1+otfVal2;
wDenom=wd;

end

