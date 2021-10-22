function Dgnd = temps_coulson(rho,GHI,surftilt,surfaz,sunel,sunaz)

    z = (90 - sunel)*pi/180;
    Dgnd = 0.5*rho*GHI.*(1-cosd(surftilt)).*(1+sin(z/2).^2).*abs(cosd(surfaz-sunaz));
end