function [ ret ] = writeApoVector( vec,otf,cutoff)
    [h,w]=size(vec);
    cnt=[h/2+1,w/2+1];
    for y=1:h
        for x=1:w
            rad=hypot(y-cnt(1),x-cnt(2));
            cycl=rad*otf.cyclesPerMicron;
            frac=cycl/(otf.cutoff*cutoff);
            if frac<0 || frac>1
                valIdealotf=0;
            else
                valIdealotf=(1/pi)*(2*acos(frac)-sin(2*acos(frac)));
            end
            vec(y,x)=valIdealotf;
        end
    end
    ret=vec;
end

