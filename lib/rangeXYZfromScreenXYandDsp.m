function [Ixyz] = rangeXYZfromScreenXYandDsp(LorRorC,IPDm,Idsp,IsmpXm,IsmpYm,screenZm)

% function [Ixyz] = rangeXYZfromScreenXYandDsp(LorRorC,IPDm,Idsp,IsmpXm,IsmpYm,screenZm)
%
%   example call: rangeXYZfromScreenXYandDsp('L',0.065,Ldsp,LsmpX,LsmpY,3.0)
%
%                 rangeXYZfromScreenXYandDsp('R',0.065,Rdsp,RsmpX,RsmpY,3.0)
%
% transforms positions of pixles in projection plane to 
% rangeXYZ by solving simultaneous parametric line equations 
% through the corresponding and nodal points
%
% LorRorC:      which image to use as reference
%            'L' -> left      image as reference
%            'R' -> right     image as reference
%            'C' -> cyclopean image as reference
% IPDm:      inter-ocular separation
% Idsp:      disparity image (pixels)                       [ n x m ]
%            (-) -> uncrossed w.r.t screen/projection-plane
%            (+) -> uncrossed w.r.t screen/projection-plane
% IsmpXm:    x-positions of pixels in screen/projection-plane (meters)
% IsmpYm:    y-positions of pixels in screen/projection-plane (meters)
%            NOTE: Idsp, IsmpXm, and IsmpYm MUST be consistent with LorRorC value
% screenZm:  z-position  of screen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ixyz:   coordinates of 3D point in cartesian coordinates assuming

if numel(screenZm)==1
   screenZm=screenZm*ones(size(IsmpXm)); 
end

% PIXEL WIDTHS IN MM
w = diff(IsmpXm(1,1:2));

if strcmp(LorRorC,'L')
    Ixyz(:,:,1) =  IPDm.*IsmpXm./(IsmpXm - (IsmpXm - IPDm - Idsp.*w));
    Ixyz(:,:,2) =  IPDm.*IsmpYm./(IsmpXm - (IsmpXm - IPDm - Idsp.*w));
    Ixyz(:,:,3) = -IPDm.*screenZm./((IsmpXm - IPDm - Idsp.*w) - IsmpXm);
elseif strcmp(LorRorC,'R')
    Ixyz(:,:,1) = -IPDm.*IsmpXm./(IsmpXm - (IsmpXm + IPDm + Idsp.*w));
    Ixyz(:,:,2) = -IPDm.*IsmpYm./(IsmpXm - (IsmpXm + IPDm + Idsp.*w));
    Ixyz(:,:,3) =  IPDm.*screenZm./((IsmpXm + IPDm + Idsp.*w) - IsmpXm);
elseif strcmp(LorRorC,'C')
    error(['rangeXYZfromScreenXYandDsp.m: WARNING! incorrect equations. FIX IT!!!']);
%     Ixyz(:,:,1) = IsmpXm;
%     Ixyz(:,:,2) = IsmpYm;
%     Ixyz(:,:,3) =  IPDm.*screenZm./((IsmpXm + IPDm + Idsp.*w) - IsmpXm);
Ixyz(:,:,3) =  IPDm.*screenZm./((IsmpXm + IPDm + Idsp.*w) - IsmpXm);
Ixyz(:,:,3) = -IPDm.*screenZm./((IsmpXm - IPDm - Idsp.*w) - IsmpXm);
else
    error(['rangeXYZfromScreenXYandDsp: WARNING! unhandled LorRorC value: ' LorRorC])
end
