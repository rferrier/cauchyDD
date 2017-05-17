function [rgb] = rgbmap (x)
% This function computes the RGB color associated to a value in [0,1]
% Values taken from GMSH

 if x<1/8
    raide = 0;
    grine = 0;
    bloue = 4*x+.5;
 elseif x<3/8
    raide = 0;
    grine = 4*x-.5;
    bloue = 1;
 elseif x<5/8
    raide = 4*x-1.5;
    grine = 1;
    bloue = -4*x+2.5;
 elseif x<7/8
    raide = 1;
    grine = -4*x+3.5;
    bloue = 0;
 else %x<1
    raide = -4*x+4.5;
    grine = 0;
    bloue = 0;
 end

 rgb = [raide,grine,bloue];
end
