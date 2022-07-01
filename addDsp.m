function [LitpRCdsp, RitpRCdsp, bIndGdFXN] = addDsp(LitpRC,RitpRC,dspDeg,LppXm,LppYm,RppXm,RppYm,IppZm,IPDm)
% XXX change IDPm to LExyz RExyz
% function [LitpRCdsp, RitpRCdsp, bIndGdFXN] =  LRSIcorrespondingPointAddDisparityExact(LitpRC,RitpRC,dspDeg,LppXm,LppYm,RppXm,RppYm,IppZm,IPDm)
%
%   example call:  % ADD UNCROSSED (-) DISPARITY TO SCENE POINT ASSOCIATED W. CORRESPONDING POINTS
%                  [LitpRCdsp, RitpRCdsp] = LRSIcorrespondingPointAddDisparityExact([144 1677.9],[144 1711.1],-15,LppXm,LppYm,RppXm,RppYm,LRSIprojPlaneDist(),LRSIcameraIPD())
%
%                  % ADD  CROSSED  (+) DISPARITY TO SCENE POINT ASSOCIATED W. CORRESPONDING POINTS
%                  [LitpRCdsp, RitpRCdsp] = LRSIcorrespondingPointAddDisparityExact([144 1677.9],[144 1711.1],+15,LppXm,LppYm,RppXm,RppYm,LRSIprojPlaneDist(),LRSIcameraIPD())
%
% WORKS VECTORIZED
% add disparity to scene point associated with corresponding image points
% while maintaining cyclopean direction
% uncrossed disparity (-) makes the target point farther than fixation
% crossed   disparity (+) makes the target point closer  than fixation
%
% LitpRC:    LE corresponding point pixel co-ordinates [ 1 x 2 ]->[ N x 2 ] (CORNER COORDINATES)
% RitpRC:    RE corresponding point pixel co-ordinates [ 1 x 2 ]->[ N x 2 ]
% dspDeg: disparity to add (in deg)
%            (-) -> uncrossed disparity
%            (+) ->   crossed disparity
% ILorR:     plane using given provided coordinate system
% IppXm:     x-positions of projection plane samples in metres in ILorR coord system
% IppYm:     y-positions of projection plane samples in metres in ILorR coord system
% IppZm:     projection plane distance in metres
% IPDm:      interpupillary distance in metres
% LitpRCdsp: LE crop-location pixel co-ordinates w. desired disparity
% RitpRCdsp: RE crop-location pixel co-ordinates w. desired disparity
% bIndGdFXN: boolean indicating if the fixation point is ahead of the eyes
%            1 -> Fixation point is in front of the eyes
%            0 -> Fixation point is behind      the eyes (BAD POINT)

    dspDeg=-dspDeg; % XXX NEED TO FIX THIS BELOW
    if abs(LitpRC(1)-RitpRC(1)) ~= 0
        warning('LRSIcorrespondingPointAddDisparityExact: Warning! LE and RE images are at different elevations. Check corresponding point');
        disp(num2str(LitpRC(1)-RitpRC(1)));
    end
    X=(1:size(LppXm,2))';
    Y=(1:size(LppYm,1))';
    LX=LppXm(1,:)';
    LY=LppYm(:,1);
    RX=RppXm(1,:)';
    RY=RppYm(:,1);
    LitpXYm(:,1) = XYZ_project.interp1qr(X,LX,LitpRC(:,2));
    LitpXYm(:,2) = XYZ_project.interp1qr(Y,LY,LitpRC(:,1));
    % RIGHT-EYE
    RitpXYm(:,1) = XYZ_project.interp1qr(X,RX,RitpRC(:,2));
    RitpXYm(:,2) = XYZ_project.interp1qr(Y,RY,RitpRC(:,1));


    % XY LOCATION IN METERS OF LE AND RE CORRESPONDING IMAGE POINTS
    % LEFT-EYE

    %LitpXYm(:,1) = interp1(1:size(LppXm,2),LppXm(1,:),LitpRC(:,2));
    %LitpXYm(:,2) = interp1(1:size(LppYm,1),LppYm(:,1),LitpRC(:,1));
    % RIGHT-EYE
    %RitpXYm(:,1) = interp1(1:size(RppXm,2),RppXm(1,:),RitpRC(:,2));
    %RitpXYm(:,2) = interp1(1:size(RppYm,1),RppYm(:,1),RitpRC(:,1));

    % ERROR CHECKING
    % if abs(LitpXYm(2)-RitpXYm(2)) ~= 0, error('LRSIcorrespondingPointAddDisparity: Warning! LE and RE images are at different elevations. Check corresponding point'); end

    % XY LOCATION IN METERS OF LE AND RE CORRESPONDING IMAGE POINTS IN CYCLOPEAN EYE COORDINATE SYSTEM
    LitpXYmC = [LitpXYm(:,1)-(IPDm/2), LitpXYm(:,2)];
    RitpXYmC = [RitpXYm(:,1)+(IPDm/2), RitpXYm(:,2)];

    % CYCLOPEAN CO-ORDINATES OF INITIAL TARGET POSITION IN 3-SPACE (USING INTERSECTION POINT)
    tgtXYZm  = intersectLinesFromPoints([ -(IPDm/2), 0 , 0 ],[ LitpXYmC(:,1), LitpXYmC(:,2), repmat(IppZm,size(LitpXYmC,1),1)],[ +(IPDm/2), 0 , 0 ],[ RitpXYmC(:,1), RitpXYmC(:,2), repmat(IppZm,size(LitpXYmC,1),1)]);
    tgtXm    = tgtXYZm(:,1);
    tgtYm    = tgtXYZm(:,2);
    tgtZm    = tgtXYZm(:,3);

    % VDisp CARTESIAN COORDINATES OF INITIAL FIXATION POSITION
    %  display(['LRSIcorrespondingPointAddDisparityExact: Fixation location in 3-space: [' num2str(tgtXYZm,'%.2f') '] in meters']);

    % CALCULATE ELEVATION ANGLE
    tgtElevDeg = atand(tgtYm./tgtZm);

    %% REPORT IF INITIAL TARGET POSITION IS 'OUT-OF-RANGE'
    bIndGdFXN = 1;  % SET bIndGdFXN = 1 BY DEFAULT

    % CHECK FOR CL BEYOND INFINITY %
    if dspDeg > 0  % ONLY A POSSIBILITY FOR CROSSED DISPARITIES
    tgtDstMtr = sqrt(sum(tgtXYZm.^2));
    maxDstMtr = (IPDm./2)/(tand((dspDeg)./2));
    % FLAG CASES WHERE FIXATION POINT IS BEHIND THE EYE
        if tgtDstMtr >= maxDstMtr
            bIndGdFXN = 0;
            % disp(['LRSIcorrespondingPointAddDisparityExact: Warning! Target distance = ', num2str(tgtDstMtr) ' m. For disparity ' num2str(dspArcMin) ' arcmin distance for parallel gaze is ' num2str(maxDstMtr) ' m. Fixation point will be behind the eyes' ]);
        end
    end

    %% TAKES INITIAL TARGET POSITION, MOVES TO ADD REQUIRED VERGENCE DEMAND AND RETURNS CYCLOPEAN CO-ORDINATES OF LE AND RE IMAGES
    % DISPARITY TO ADD IN DEGREE

    % VERGENCE OF INITIAL TARGET POSITION (IN EPIPOLAR PLANE)
    vrgCPdeg   = acosd((tgtXm.^2 + tgtYm.^2 + tgtZm.^2 - (IPDm./2).^2)./((sqrt((tgtXm + (IPDm./2)).^2 + tgtYm.^2 + tgtZm.^2)).*(sqrt((tgtXm - (IPDm./2)).^2 + tgtYm.^2 + tgtZm.^2))));

    % VERSION ANGLE OF INITIAL TARGET POSITION (IN EPIPOLAR PLANE)
    vrsCPdeg = -sign(tgtXm).*acosd(sqrt((tgtYm.^2 + tgtZm.^2)./(tgtXm.^2 + tgtYm.^2 + tgtZm.^2)));

    % TOLERANCE FOR CHECKING IF THE VERSION ANGLE  IS CLOSE TO ZERO
    tolVrsDeg = 0.001;
    tolVrgDeg = 1/60;

    % CYLOPEAN  XY CO-ORDINATES OF MOVED TARGET IMAGE (IN THE ORIGINAL SPACE)

    if abs(vrgCPdeg - dspDeg) < tolVrgDeg                                                                 % DEALS WITH CASE OF PARALLEL GAZE
        % IF GAZE IS PARALLEL, IMAGE PLANE CO-ORDINATES ARE OBTAINED WITHOUT SOLVING FOR THE NEW FIXATION-POINT CO-ORDINATES IN 3-SPACE
        LitpXYmCdsp = [ -(IPDm./2) - (IppZm./cosd(tgtElevDeg)).*tand(vrsCPdeg), LitpXYmC(:,2)];
        RitpXYmCdsp = [ +(IPDm./2) - (IppZm./cosd(tgtElevDeg)).*tand(vrsCPdeg), RitpXYmC(:,2)];
    else % DEALS WITH CASE WHERE FIXATION LOCATION IS AT FINITE RANGE
        if abs(vrsCPdeg)<tolVrsDeg % ZERO VERSION ANGLE IN EPIPOLAR PLANE
            XclTrnsMtr = 0;
            ZclTrnsMtr = (IPDm./2)./(tand((vrgCPdeg - dspDeg)./2));
        else    % NONZERO VERSION ANGLE IN EPIPOLAR PLANE
            if bIndGdFXN == 1
                % TARGET AHEAD OF THE EYES IN THE EPIPOLAR PLANE
                ZclTrnsMtr = ((cosd(vrsCPdeg)).^2).*(IPDm./2).*( cotd(vrgCPdeg - dspDeg) + sqrt( (cotd(vrgCPdeg - dspDeg)).^2+(secd(vrsCPdeg)).^2 ) ); % POSITIVE ROOT IS CHOSEN OF THE QUADRATIC IN ZclTrnsMtr
            elseif bIndGdFXN == 0
                % TARGET BEHIND THE EYES IN THE EPIPOLAR PLANE
                ZclTrnsMtr = ((cosd(vrsCPdeg)).^2).*(IPDm./2)*( cotd(vrgCPdeg - dspDeg) - sqrt( (cotd(vrgCPdeg - dspDeg)).^2+(secd(vrsCPdeg)).^2 ) ); % NEGATIVE ROOT IS CHOSEN OF THE QUADRATIC IN ZclTrnsMtr
            end
            XclTrnsMtr = -tand(vrsCPdeg).*ZclTrnsMtr;
        end
        %IMAGE PLANE INTERSECTION XY CO-ORDINATES FOR TARGET IN FRONT
        LitpXYmCdsp = [ -(IPDm./2) + (IppZm./(ZclTrnsMtr.*cosd(tgtElevDeg))).*(XclTrnsMtr + (IPDm./2))  , LitpXYmC(:,2)];
        RitpXYmCdsp = [ +(IPDm./2) + (IppZm./(ZclTrnsMtr.*cosd(tgtElevDeg))).*(XclTrnsMtr - (IPDm./2))  , RitpXYmC(:,2)];
    end

    %% TAKES CYCLOPEAN XY IMAGE CO-ORDINATES AND RETURNS LE AND RE IMAGE CO-ORDINATES FOR THE REQUIRED EYE
    LitpXYmDsp = LitpXYmCdsp + [+(IPDm/2),0];
    RitpXYmDsp = RitpXYmCdsp + [-(IPDm/2),0];

    %% TAKES LE AND RE IMAGE PLANE CO-ORDINATES IN METERS AND RETURNS PIXEL CO-ORDINATES
    %LitpRCdsp(:,1) = interp1(LppYm(:,1),1:size(LppYm,1),LitpXYmDsp(:,2));
    %LitpRCdsp(:,2) = interp1(LppXm(1,:),1:size(LppXm,2),LitpXYmDsp(:,1));
    %RitpRCdsp(:,1) = interp1(RppYm(:,1),1:size(RppYm,1),RitpXYmDsp(:,2));
    %RitpRCdsp(:,2) = interp1(RppXm(1,:),1:size(RppXm,2),RitpXYmDsp(:,1));
    LitpRCdsp(:,1) = XYZ_project.interp1qr(flip(LY,1),flip(Y,1),LitpXYmDsp(:,2));
    LitpRCdsp(:,2) = XYZ_project.interp1qr(LX,X,LitpXYmDsp(:,1));
    RitpRCdsp(:,1) = XYZ_project.interp1qr(flip(RY,1),flip(Y,1),RitpXYmDsp(:,2));
    RitpRCdsp(:,2) = XYZ_project.interp1qr(RX,Y,RitpXYmDsp(:,1));
end
