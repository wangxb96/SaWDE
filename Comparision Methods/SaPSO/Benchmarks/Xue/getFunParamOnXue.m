function [  XRmin, XRmax, Lbound, Ubound,VTR, maxfeval ] = getFunParamOnWang( D, fun )
    XRmin =0*ones(1,D);
    XRmax =1*ones(1,D);
    Lbound = XRmin;
    Ubound = XRmax;
    VTR = 0;
    maxfeval=1000000;




