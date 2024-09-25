function [ret] = applyOtf( vec,otf,band,kx,ky,useAttenuation,write )
    ret=otfToVector(vec,otf,band,kx,ky,useAttenuation,write);
end

