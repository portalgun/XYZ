function [AitpRC,BitpRC,CitpRC,CsmpErrDeg,DvrgDffDeg,indNaN] = LRSIcorrespondingPointA2B(LorR,ActrRC,Axyz,Bxyz,AppXm,AppYm,bPLOT,bNaNhide)
% function [LitpRC,RitpRC,CitpRC,CsmpErrDeg,DvrgDffDeg,indNaN] = LRSIcorrespondingPointA2B(ActrRC,Lxyz,Rxyz,LppXm,LppYm,bPLOT,bNaNhide)
%
%   example call:   [~,~,~,~,Lxyz,Rxyz,LppXm,LppYm,RppXm,RppYm] = loadLRSIimage(62,1,1,'PHT','img');
%                   [AitpRC,BitpRC,CitpRC,CsmpErrDeg,DvrgDffDeg,indNaN] = LRSIcorrespondingPointA2B('L',[628 1048; 420 360],Lxyz,Rxyz,LppXm,LppYm,0,0);
%
% from point in eye A image, finds corresponding points in eye B image. eg. if A == L, then B == R
%
% ActrRC:        coordinates of points in anchor eye for which to
%                find other eye corresponding points
% Lxyz:          left  eye cartesian coordinates to points in scene
%                (in meters)
% Rxyz:          right eye cartesian coordinates to points in scene
%                (in meters)
% AppXm:         x-positions in anchor eye of pixel center in projection plane
% AppYm:         y-positions in anchor eye of pixel center in projection plane
% bPLOT:         1 -> plot
%                0 -> not
% bNaNhide:      Keeps values that should be NaNs as one, for further
%                computation on entire matrix outside of the
%                function.  Use indNaN as reference to change to NaN
%                later.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AitpRC:        interpolated corresponding points in anchor eye's image
% RitpRC:        interpoalted corresponding points in other  eye's image
% CitpRC:        Interpolated corresponding points in cyclopean eye's image
% CsmpErrDeg:    vergence error of 'sampled' cyclopean 3D sampled
%                points assessed by examining vergence demand
%                differences between 'sampled' and 'interpolated'
%                cyclopean points.  this error should not exceed the
%                angular spacing of the samples in the projection plane
% DvrgDffDeg:    vergence difference between sampled eye A and eye B points
% indNaN:        index of corresponding points that are NaN
% --------------------------------------------------------------------
%TODO
%
%vergenceFromRangeXYZ - ok
%intersectLinePlane   - ok
%intesecitonPoint     - ok
%createLine3d         - ok
%combine L2R and R2L  - ok
%Plot portion         - test

% -------------------------------------------------------------------------------
% CHECKS
%test mode returns testing variables to workspace for comparison with non-vectorized version of code
if ~exist('bPLOT','var') || isempty(bPLOT); bPLOT = 0; end
if ~exist('bNaNhide','var') || isempty(bNaNhide); bNaNhide=0; end

if any(mod(ActrRC(:,1),1)) && strcmp(LorR,'L')
    error(['LRSIcorrespondingPointL2R: WARNING! ' ActrRC(:,1) ' must be integer valued']);
end
% -------------------------------------------------------------------------------
% CONSTANTS
AppZm           = LRSIprojPlaneDist(1);
IPDm            = LRSIcameraIPD(1);

IszRC=[size(Axyz,1) size(Axyz,2)];
res=size(Axyz,1)*size(Axyz,2)*size(Axyz,3);

N=size(ActrRC,1);
[n,m,~]=size(Axyz);
PrjPln  = createPlane([0 0 AppZm], [1 0 AppZm], [0 1 AppZm]);
distYperPix=(AppYm(end,end)-AppYm(1,1))/size(AppYm,1);
distXperPix=(AppXm(end,end)-AppXm(1,1))/size(AppXm,2);
pixPerMtr    = size(AppXm,2)./diff([AppXm(1) AppXm(end)]);

if strcmp(LorR,'L') %Assign L and R to a and B according to starting pont
    A='L';
    B='R';
else
    A='R';
    B='L';
end


% XYZ COORDINATES OF PROJ PLANE & LEFT/RIGHT/CYCLOPEAN EYE (in left eye coordinate frame)
%L
if A=='L'
    AxyzEye     = [    0   0 0];
    BxyzEye     = [+IPDm   0 0];
    CxyzEye     = [+IPDm/2 0 0];
    % XminXmax    = minmax([Lxyz(:,:,1) Rxyz(:,:,1)]); % for robustness (see below)
    % ZminZmax    = minmax([Lxyz(:,:,3) Rxyz(:,:,3)]); % for robustness (see below)
elseif A=='R'
    AxyzEye     = [    0   0 0];
    BxyzEye     = [-IPDm   0 0];
    CxyzEye     = [-IPDm/2 0 0];
    % XminXmax    = minmax([Lxyz(:,:,1) Rxyz(:,:,1)]); % for robustness (see below)
    % ZminZmax    = minmax([Lxyz(:,:,3) Rxyz(:,:,3)]); % for robustness (see below)
end

bDebug=false;
if bDebug
    disp('--------------------------')
    [~,~,Xctr]=Map.getCtrVal(AppXm);
    [~,~,Yctr]=Map.getCtrVal(AppYm);
    %[~,~,Mctr]=Map.getCtrVal(M);
    [AU,AL,Actr]=Map.getCtrVal(Axyz);

    disp([ ...
            'LorR    ' LorR newline ...
            'BxyzEye ' Num.toStr(BxyzEye) newline ...
            'Xctr    ' num2str(Xctr) newline ...
            'AU      ' num2str(AU) newline ...
            'AL      ' num2str(AL) newline ...
            'Actr    ' num2str(Actr) newline ...
            'dX      ' num2str(distXperPix) newline ...
            'dY      ' num2str(distYperPix) newline ...
            'IPDm    ' num2str(IPDm) newline ...
            'Z       ' num2str(AppZm) newline ...
            ]);
end


% -------------------------------------------------------------------------------
%FIND NaNs in ActrRC
[indNaN,~] = find(isnan(ActrRC) | ActrRC<=0 | ActrRC > repmat(IszRC,size(ActrRC,1),1)); %Handle NaNs in ActrRC
if ~isempty(indNaN) %temporarily assign 1s to index values, so that the entire matrix computes
    ActrRC(indNaN,:)=1;
end
% -------------------------------------------------------------------------------


% XYZ AT CENTER OF LE SAMPLED PATCHES IN 3D
try
    ind=sub2ind(size(Axyz),ActrRC(:,1),ActrRC(:,2));
catch ME
    size(Axyz)
    ActrRC
    rethrow(ME)
end
ind=[ind ind+n*m ind+2*n*m ];

%POINTS OUT OF RANGE
invInd=any(ind>res,2);
ind(invInd,:)=1; %temorparily assign 1

AxyzCtr = Axyz(ind);

%% Handle NaNs in AxyzCtr
[ind,~]    = find(isnan(AxyzCtr));
indNaN     = [ind; indNaN];
[ind,~]    = find(invInd);
indNaN     = [ind; indNaN];
indNaN=unique(indNaN);
if ~isempty(indNaN) %temporarilty assign 1s to index values, so that the entire matrix computes
    AxyzCtr(indNaN,:)=1;
    N=size(ActrRC,1);
end

%Remove NaN's

%%% CHECKED
% VERGENCE ANGLE AT SAMPLED SCENE POINT
AVrgCtrDeg = vergenceFromRangeXYZVec(A,IPDm,AxyzCtr);

% LE AND RE LINES OF SIGHT TO LEFT EYE POINT
A2AvctLOS = createLine3d(AxyzEye,AxyzCtr);
A2BvctLOS = createLine3d(BxyzEye,AxyzCtr);

% INTERSECTION OF PROJECTION PLANE (IPP) W. LINES OF SIGHT TO POINT IN L IMAGE
A2AxyzIPP = intersectLinePlane(A2AvctLOS,PrjPln); % must be    center of pixel location
A2BxyzIPP = intersectLinePlane(A2BvctLOS,PrjPln); % may not be center of pixel location


% -------------------------------------------------------------------------------
% 3D SAMPLED POINT IN IMAGE B NEAREST THE CORRESPONDING POINT (INDEX,ROW,COL)
BctrRC=zeros(N,2);
BctrRC0=zeros(N,2);

%FASTEST
closestY = interp1(AppYm(:,1),AppYm(:,1),A2BxyzIPP(:,2),'nearest','extrap');
closestX = interp1(AppXm(1,:),AppXm(1,:),A2BxyzIPP(:,1),'nearest','extrap');

BctrRC(:,1)=ceil((closestY./distYperPix)+IszRC(1)/2);

BctrRC(:,2)=ceil((closestX-CxyzEye(1))./distXperPix+IszRC(2)/2);

ind=(BctrRC(:,1)==0 | BctrRC(:,2)==0);
BctrRC(ind,1)=1;
BctrRC(ind,2)=1;
%%%


% CORRESPONDING POINTS MUST HAVE SAME VERTICAL VALUE (i.e. THEY LIE IN AN EPIPOLAR PLANE)
ind=find(ActrRC(:,1) ~= BctrRC(:,1));
if ~isempty(ind)
    disp(['LRSIcorrespondingPointL2R: WARNING! forcing vertical pixel row to be the same.']);
    BctrRC(ind,:)=[ActrRC(ind,1) BctrRC(ind,2)];
end

%XXX
%FIND CPs OUTSIDE OF RANGE
badInd=find(BctrRC(:,2)>IszRC(2) | BctrRC(:,2)<=0);

%%% CHECKED BELOW
BctrRC(badInd,2)=1;
indNaN=unique([indNaN; badInd]);

% XYZ OF RE SAMPLED POINT NEAREST THE TRUE CORRESPONDING POINT IN PROJECTION PLANE (in left coordinate frame)
M=bsxfun(@plus,Bxyz,reshape(BxyzEye,1,1,3));
BxyzCtr=zeros(N,3);
sub=[repelem(BctrRC(:,1),3,1), repelem(BctrRC(:,2),3,1), repmat([1:3]',size(BctrRC,1),1)];
ind=sub2ind(size(M),sub(:,1),sub(:,2),sub(:,3));

BxyzCtr = M(ind);
BxyzCtr = reshape(BxyzCtr,3,size(BxyzCtr,1)/3)';

%OLD
%M=bsxfun(@plus,permute(Bxyz,[2,3,1]),BxyzEye);
%BBxyzCtr=zeros(N,3);
%for i = 1:N
%    BBxyzCtr(i,:) = M(BctrRC(i,2),:,BctrRC(i,1)); % + RxyzEye puts BxyzCtr in LE coordinate system
%end
%[BBxyzCtr(n,:)  BxyzCtr(n,:)]


%%% CHECKED BELOW
% VERGENCE ANGLE AT SAMPLED SCENE POINT FOR OPPOSITE EYE
BVrgCtrDeg   = vergenceFromRangeXYZVec(A,IPDm,reshape(BxyzCtr,[N 1 3]));

% LE AND RE LINES OF SIGHT TO RIGHT EYE POINT
% R2LvctLOS    = createLine3d(LxyzEye,BxyzCtr);
% R2RvctLOS    = createLine3d(RxyzEye,BxyzCtr);

% INTERSECTION OF PROJECTION PLANE (IPP) W. LINES OF SIGHT TO POINT IN R IMAGE
B2AxyzIPP    = intersectLinePlane(createLine3d(AxyzEye,BxyzCtr),PrjPln); % must be    center of pixel location
B2BxyzIPP    = intersectLinePlane(createLine3d(BxyzEye,BxyzCtr),PrjPln); % may not be center of pixel location

% SAMPLED 3D POINT NEAREST THE CORRESPONDING POINT IN 3D (may not correspond to 3D surface)
CxyzCtr      = intersectLinesFromPoints(AxyzEye,AxyzCtr,BxyzEye,BxyzCtr,1000);      % 3D coordinates  of   sampled    point
%CxyzCtr      = intersectionPointVec(AxyzEye,AxyzCtr,BxyzEye,BxyzCtr);      % 3D coordinates  of   sampled    point
CvrgCtrDeg   = vergenceFromRangeXYZVec(A,IPDm,reshape(CxyzCtr,[N 1 3])); % vergence demand of   sampled    point

% INTERPOLATED 'TRUE' POINT (if LE and RE sample are on surface, ITP point should be very near to surface)
CxyzItp      = intersectLinesFromPoints(CxyzEye,CxyzCtr,AxyzCtr,BxyzCtr,1000);      % 3D coordinates  of interpolated point
%CxyzItp      = intersectionPointVec(CxyzEye,CxyzCtr,AxyzCtr,BxyzCtr);      % 3D coordinates  of interpolated point

% ERROR CHECKING... INTERPOLATED POINT CANNOT BE NEARER/FARTHER THAN MIN/MAX DISTANCE
% if CxyzItp(1) < XminXmax(1) || CxyzItp(1) > XminXmax(2) || CxyzItp(3) < ZminZmax(1) || CxyzItp(3) > ZminZmax(2),
%     CxyzItp = BxyzCtr; disp(['LRSIcorrespondingPointL2R: WARNING! forcing CxyzItp = BxyzCtr']);
% end

%% VERGENCE ANGLE
CvrgItpDeg   = vergenceFromRangeXYZVec(A,IPDm,reshape(CxyzItp,[N 1 3])); % vergence demand of interpolated point


% INTERSECTION OF PROJECTION PLANE (IPP) W. LINES OF SIGHT FROM LE, RE, CE TO INTERPOLATED CYCLOPEAN POINT
C2AxyzIPP    = intersectLinePlane(createLine3d(AxyzEye,CxyzItp),PrjPln);
C2BxyzIPP    = intersectLinePlane(createLine3d(BxyzEye,CxyzItp),PrjPln);
C2CxyzIPP    = intersectLinePlane(createLine3d(CxyzEye,CxyzItp),PrjPln);


% LE AND RE IMAGE SHIFTS IN METERS FOR 'TRUE' INTERPOLATED (ITP) POINT TO NULL ERROR IN 3D SAMPLED POINT
AitpShftXm   = C2AxyzIPP(:,1) - A2AxyzIPP(:,1);
BitpShftXm   = C2BxyzIPP(:,1) - B2BxyzIPP(:,1);
CitpShftXm   = C2CxyzIPP(:,1) - A2AxyzIPP(:,1);


% LE AND RE IMAGE SHIFTS IN PIXELS FOR 'TRUE' INTERPOLATED (ITP) POINT TO NULL ERROR IN 3D SAMPLED POINT
AitpShftXpix = AitpShftXm.*pixPerMtr;
BitpShftXpix = BitpShftXm.*pixPerMtr;
CitpShftXpix = CitpShftXm.*pixPerMtr;


% INTERPOLATED LOCATION OF CORRESPONDING POINT (IN PIXELS)
AitpRC    = ActrRC + [zeros(N,1) AitpShftXpix];
BitpRC    = BctrRC + [zeros(N,1) BitpShftXpix];
CitpRC    = ActrRC + [zeros(N,1) CitpShftXpix];


%% DISPARITY BETWEEN BEST 3D SAMPLED POINT AND BEST INTERPOLATED 3D POINT
CsmpErrDeg = ( CvrgItpDeg - CvrgCtrDeg );

% SAMPLED FIXATION DISTANCE (from cyclopean eye
CdstErrM   = sqrt(sum(abs((bsxfun(@minus,CxyzItp,CxyzEye)).^2),2));
AV30=CdstErrM;

% VERGENCE DIFFERENCE BETWEEN NEAREST SAMPLED POINTS (SHOULD BE TINY)
DvrgDffDeg = AVrgCtrDeg - BVrgCtrDeg;
AV31=DvrgDffDeg;

% HANDLE NaN CORRESPONDING POINTS -> 1, changed to NaN in parent function
if ~isempty(indNaN)
    if bNaNhide == 1
        AitpRC(indNaN,:)=1;
        BitpRC(indNaN,:)=1;
        CitpRC(indNaN,:)=1;
    else
        AitpRC(indNaN,:)=NaN;
        BitpRC(indNaN,:)=NaN;
        CitpRC(indNaN,:)=NaN;
    end
    CsmpErrDeg(indNaN,:) = 100; DvrgDffDeg(indNaN,:) = 100;
end

% -------------------------------------------------------------------------------
% PLOT
if bPLOT == 1
    %%
    figure(11111);
    set(gcf,'position',[120    15   615   627]);
    hold on;
    % axis([max([abs(minmax([2*[-IPDm + IPDm] + minmax(LxyzCrp(:,:,1))]))  abs(minmax([0 + minmax(LxyzCrp(:,:,2))]))])*[-1 1 -1 1] minmax([0 minmax(LxyzCrp(:,:,3))])])
    % axis([-.1 .1 -.1 .1 0  10.5])

    if A=='L'
        plot3([0 IPDm],[0 0],[ 0 0],'ko');
    else
        plot3([AxyzEye(1) BxyzEye(1)],[AxyzEye(2) BxyzEye(2)],[AxyzEye(3) BxyzEye(3)],'ko');
    end
    box on; grid on;
    Fig.format('X','Y');
    zlabel('Z');

    % LEFT IS RED, RIGHT IS BLUE
    if A=='L'
        a='r';
        b='b';
    else
        b='r';
        a='b';
    end
    ao=[a 'o'];
    bo=[b 'o'];

    %PLOT
    for i = 1:N
        AR=ActrRC(i,1)+[-5:5];
        AC=ActrRC(i,2)+[-5:5];
        BR=BctrRC(i,1)+[-5:5];
        BC=BctrRC(i,2)+[-5:5];
        if any(AR < 1) || any(AC < 1) || any(BR < 1) || any(BC < 1) || any(AR > IszRC(1)) || any(AC > IszRC(2)) || any(BR > IszRC(1))  || any(BC > IszRC(2))
            continue
        end
        plot3(AxyzCtr(i,1),AxyzCtr(i,2),AxyzCtr(i,3),ao,'markersize',12,'linewidth',2);
        plot3(BxyzCtr(i,1),BxyzCtr(i,2),BxyzCtr(i,3),bo,'markersize',12,'linewidth',2);
        plot3(CxyzCtr(i,1),CxyzCtr(i,2),CxyzCtr(i,3),'ko','markersize',12,'linewidth',2);
        plot3(CxyzItp(i,1),CxyzItp(i,2),CxyzItp(i,3),'mo','markersize',12,'linewidth',2);
        plot3(...
                Axyz(AR,AC,1),...
                Axyz(AR,AC,2),...
                Axyz(AR,AC,3),a);
        plot3(...
                Bxyz(BR,BC,1)+IPDm,...
                Bxyz(BR,BC,2),...
                Bxyz(BR,BC,3),b);

        %drawLine3d(createLine3d(BxyzCtr(i,:),AxyzEye),b);
        %drawLine3d(createLine3d(BxyzCtr(i,:),BxyzEye),b);
        %drawLine3d(createLine3d(AxyzCtr(i,:),AxyzEye),a);
        %drawLine3d(createLine3d(AxyzCtr(i,:),BxyzEye),a);
        %drawLine3d(createLine3d(CxyzEye(i,:),CxyzCtr(i,:),'k'));
        %drawLine3d(createLine3d(AxyzCtr(i,:),BxyzCtr(i,:),'m'));

        % PLOT INTERSECTION WITH PROJECTION PLANE %XXX
        plot3(  A2AxyzIPP(i,1),  A2AxyzIPP(i,2),  A2AxyzIPP(i,3),ao,'markersize',12,'linewidth',2);
        plot3(  A2BxyzIPP(i,1),  A2BxyzIPP(i,2),  A2BxyzIPP(i,3),ao,'markersize',12,'linewidth',2);
        plot3(C2AxyzIPP(i,1),C2AxyzIPP(i,2),C2AxyzIPP(i,3),'mo','markersize',12,'linewidth',2);
        plot3(C2BxyzIPP(i,1),C2BxyzIPP(i,2),C2BxyzIPP(i,3),'mo','markersize',12,'linewidth',2);
        %drawLine3d(createLine3d(AxyzEye,CxyzItp),'m');
        %drawLine3d(createLine3d(BxyzEye,CxyzItp),'m');
    end

    % PROJECTION PLANE LINE AT EYE HEIGHT
    drawLine3d(createLine3d([-1 0 3],[1 0 3]),'k');
    % axis([0.0075 .0275 -.01 .01 3.44  3.5])
    set(gca,'ydir','reverse')
    view([0 0])
    killer = 1;
end
