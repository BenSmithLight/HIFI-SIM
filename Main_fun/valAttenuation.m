function [va]= valAttenuation(dist,str,fwhm)
    va=1-str*(exp(-power(dist,2)/(power(0.5*fwhm,2))));
end

