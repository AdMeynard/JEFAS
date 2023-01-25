function [psih, dpsih, ddpsih] = cmor(Fb,Fc,xi)
    psih = sqrt(Fb) * exp(-Fb^2*pi*(xi-Fc).^2);
    dpsih = -2*Fb^2*pi*(xi-Fc).*psih;
    %ddpsih = (-2*Fb^2*pi + (-2*Fb^2*pi*(xi-Fc)).^2).*psih;
    ddpsih = 2*pi*Fb^2*(2*pi*Fb^2*Fc^2-4*pi*Fb^2*Fc.*xi+2*pi*Fb^2*xi.^2-1).*psih;
end