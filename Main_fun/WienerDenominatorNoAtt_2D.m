function [ wDenom ] = WienerDenominatorNoAtt_2D(sp,d,b)
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
        otfVal1=otfVal1/2;
        otfVal2=otfVal2/2;
    end
        wDenom=otfVal1+otfVal2;
end

