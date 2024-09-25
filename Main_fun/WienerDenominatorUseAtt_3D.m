function [ wDenom ] = WienerDenominatorUseAtt_3D(sp,d,b)

w=2*sp.imgSize;
h=2*sp.imgSize;
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
    otfVal1=otfVal1.*valAttenuation(rad1,1*sp.OtfProvider.attStrength/1.05,1*sp.OtfProvider.attFWHM)/2;
    otfVal2=otfVal2.*valAttenuation(rad2,1*sp.OtfProvider.attStrength/1.05,1*sp.OtfProvider.attFWHM)/2;
elseif b==2
%     otfVal1=otfVal1.*valAttenuation(rad1,sp.OtfProvider.attStrength/1.05,sp.OtfProvider.attFWHM)/1.0;
%     otfVal2=otfVal2.*valAttenuation(rad2,sp.OtfProvider.attStrength/1.05,sp.OtfProvider.attFWHM)/1.0;
    otfVal1=otfVal1.*valAttenuation(rad1,sp.OtfProvider.attStrength/1.00,sp.OtfProvider.attFWHM)/1.0;
    otfVal2=otfVal2.*valAttenuation(rad2,sp.OtfProvider.attStrength/1.00,sp.OtfProvider.attFWHM)/1.0;
    
else
%     otfVal1=otfVal1.*valAttenuation(rad1,sp.OtfProvider.attStrength/1.05,sp.OtfProvider.attFWHM)/1.0;
%     otfVal2=otfVal2.*valAttenuation(rad2,sp.OtfProvider.attStrength/1.05,sp.OtfProvider.attFWHM)/1.0;
    otfVal1=otfVal1.*valAttenuation(rad1,sp.OtfProvider.attStrength/1.0,sp.OtfProvider.attFWHM)/1.0;
    otfVal2=otfVal2.*valAttenuation(rad2,sp.OtfProvider.attStrength/1.0,sp.OtfProvider.attFWHM)/1.0;
end
    wDenom=otfVal1+otfVal2;
end

