classdef CPs_gen < handle
methods(Static)

    function [outL,outR]=genImgBoth(database,I,ActrRC,bAndChk,db)
    %function [outL,outR]=CPs.genImgBoth(database,I, OPTIONAL: ActrRC,bAndChk, SPEED_OPT: db);
        if nargin < 6
            db=dbInfo(database);
            if nargin < 5
                bAndChk=true;
                if nargin < 4
                    ActrRC=[];
                end
            end
        end
        if isempty(ActrRC)
            ActrRC=CPs.getAllRC(db.IszRC);
        end

        [OptsL]=CPs_gen.get_db_params(database,'L');
        [OptsR]=CPs_gen.get_db_params(database,'R');
        [Lxyz,Lm,LBxyz] = CPs_gen.get_img_params(I,OptsL);
        [Rxyz,Rm,RBxyz] = CPs_gen.get_img_params(I,OptsR);

        %isequal(OptsL.AxyzEye,OptsR.BxyzEye)
        %isequal(OptsL.BxyzEye,OptsR.AxyzEye)
        %isequal(OptsL.BppXm,OptsR.AppXm)
        %isequal(OptsL.AppXm,OptsR.BppXm)
        %isequal(OptsL.AppYm,OptsR.AppYm)
        %~isequal(Lxyz,Rxyz)
        %~isequal(Lm,Rm)

        if bAndChk
            outL=CPs_gen.get_AB_and_chk(ActrRC,Lxyz,Rxyz,Lm,Rm,OptsL,OptsR);
            outR=CPs_gen.get_AB_and_chk(ActrRC,Rxyz,Lxyz,Rm,Lm,OptsR,OptsL);
        else
            outL=struct();
            outR=struct();
            [outL.AitpRC, outL.BitpRC, outL.CitpRC, ...
             outL.indNan...
             outL.LitpRC, outL.RitpRC, ...
             outL.CsmpErrDeg, outL.dvrgDffDeg ] = CPs_gen.get_AB(ActrRC,Lxyz,Lm,OptsL,false);

            [outR.AitpRC, outR.BitpRC, outR.CitpRC, ...
             outL.indNan...
             outR.LitpRC, outR.RitpRC, ...
             outR.CsmpErrDeg, outR.dvrgDffDeg ] = CPs_gen.get_AB(ActrRC,Rxyz,Rm,OptsR,false);
        end
    end
end
methods(Static,Hidden)
    function M=get_M(BxyzEye,Bxyz)
        M=bsxfun(@plus,Bxyz,reshape(BxyzEye,1,1,3));
        %BxyzEye2(1,1,:)=BxyzEye;
        %M=Bxyz+BxyzEye2;
    end
    function [Axyz,Am,Bxyz]=get_img_params(I,Opts)
        Axyz=dbImg.getImg(Opts.database,'img','xyz',I,Opts.LorR);
        Bxyz=dbImg.getImg(Opts.database,'img','xyz',I,Opts.nLorR);
        Am=CPs_gen.get_M(Opts.BxyzEye,Bxyz);
    end
    function [Opts]=get_db_params(database,LorR, db)
        if nargin < 4
            db=dbInfo(database);
        end

        Opts=struct();
        Opts.database=db.database;
        Opts.LorR   =LorR;
        [Opts.k,Opts.nk]=CPs.getK(LorR);
        Opts.nLorR  =CPs.NOTLANDR(Opts.nk);
        Opts.IPDm   =db.IPDm;
        Opts.AxyzEye=db.(LorR).AExyz;
        Opts.BxyzEye=db.(LorR).BExyz;
        Opts.CxyzEye=db.(LorR).CExyz;
        Opts.CppXm  =db.IppXm{3};
        if LorR == 'L'
            Opts.AppXm  =db.IppXm{1};
            Opts.BppXm  =db.IppXm{2};

            Opts.AppYm  =db.IppYm{1};
        elseif LorR =='R'
            Opts.AppXm  =db.IppXm{2};
            Opts.BppXm  =db.IppXm{1};

            Opts.AppYm  =db.IppYm{2};
        end
        Opts.AppZm=db.IppZm;
    end
end
methods(Static,Access={?CPs,?CPsDB})
    function out=get_AB_and_chk(ActrRC,Axyz,Bxyz,AM,BM,OptsA,OptsB)

        [out.AitpRC,  out.BitpRC, out.CitpRC,...
         out.indNan, ...
         out.LitpRC,  out.RitpRC, ...
         out.CsmpErrDeg, out.dvrgDffDeg]                  = CPs_gen.get_AB(           ActrRC,Axyz,Bxyz,AM,OptsA,true);

        [out.BitpRCchk, out.AitpRCchk, out.CitpRCchk, ...
         out.indNanChk, ...
         out.LitpRCchk,  out.RitpRCchk]                  = CPs_gen.get_AB(round(out.BitpRC),Bxyz,Axyz,BM,OptsB,true);


        out.k=OptsA.LorR;
        out.ActrRC=ActrRC;
        out.AitpRC(out.indNan,:)=nan;
        out.BitpRC(out.indNan,:)=nan;
        out.CitpRC(out.indNan,:)=nan;
        out.LitpRC(out.indNan,:)=nan;
        out.RitpRC(out.indNan,:)=nan;

        out.AitpRCchk(out.indNanChk,:)=nan;
        out.BitpRCchk(out.indNanChk,:)=nan;
        out.CitpRCchk(out.indNanChk,:)=nan;
        out.LitpRCchk(out.indNanChk,:)=nan;
        out.RitpRCchk(out.indNanChk,:)=nan;

        out.AHOccInd=any(abs(out.AitpRC-out.AitpRCchk)>.1,2);
        out.BHOccInd=any(abs(out.BitpRC-out.BitpRCchk)>.1,2);
        %out.GdInd=out.dvrgDffDeg < .1 & out.CsmpErrDeg  < .1 & ~any(isnan(out.AitpRC),2);

        %sum(out.GdInd)
        %numel(out.GdInd)

    end
    function [AitpRC, BitpRC, CitpRC, indNan, LitpRC, RitpRC, CsmpErrDeg, DvrgDffDeg] = get_AB(ActrRC,Axyz,Bxyz,AM,Opts,bNanHide)
        bLRSI=false; %NOTE
        if ~bLRSI
            [AitpRC, BitpRC, CitpRC, indNan, CsmpErrDeg, DvrgDffDeg] = CPs_gen.get_AB_main(Opts.LorR, ActrRC,Axyz, Bxyz, ...
                                                                                   Opts.AppXm, Opts.AppYm, Opts.AppZm, ...
                                                                                   Opts.IPDm, Opts.AxyzEye, Opts.BxyzEye, Opts.CxyzEye, bNanHide);
        else
            [AitpRC,BitpRC,CitpRC,CsmpErrDeg,DvrgDffDeg,indNan] = LRSIcorrespondingPointA2B(Opts.LorR,ActrRC,Axyz,Bxyz,Opts.AppXm,Opts.AppYm,0,1);
        end
        if Opts.LorR=='L'
            LitpRC=AitpRC;
            RitpRC=BitpRC;
        elseif Opts.LorR=='R'
            LitpRC=BitpRC;
            RitpRC=AitpRC;
        end

    end
end
methods(Static,Access=private)
    function [AitpRC,BitpRC,CitpRC,indNaN,CsmpErrDeg,DvrgDffDeg] = get_AB_main(LorR, ActrRC,Axyz,Bxyz, AppXm,AppYm,AppZm,IPDm,AxyzEye,BxyzEye,CxyzEye ,bNaNhide, bDebug)

        % CONSTANTS
        IszRCD=size(Axyz);
        IszRC=IszRCD(1:2);
        res=prod(IszRCD);

        N=size(ActrRC,1);
        [n,m,~]=size(Axyz);
        PrjPln  = createPlane([0 0 AppZm], [1 0 AppZm], [0 1 AppZm]);
        distYperPix=(AppYm(end,end)-AppYm(1,1))/size(AppYm,1);
        distXperPix=(AppXm(end,end)-AppXm(1,1))/size(AppXm,2);
        pixPerMtr    = size(AppXm,2)./diff([AppXm(1) AppXm(end)]);

        %% DEBUG
        if nargin < 13
            bDebug=false;
        end
        bDebug=false; % NOTE
        if bDebug
            disp('--------------------------')
            [~,~,Xctr]=Map.getCtrVal(AppXm);
            [~,~,Yctr]=Map.getCtrVal(AppYm);
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
                    %[~,~,Mctr]=Map.getCtrVal(M);
                   %'Mctr    ' num2str(Mctr) newline ...
                   %newline ...
                   %'AxyzEye ' Num.toStr(AxyzEye) newline ...
                   %'CxyzEye ' Num.toStr(CxyzEye) newline ...
                   %'Yctr    ' num2str(Yctr) newline ...
        end

        A=LorR;

%% HANDLE NaNs
        indNaN=CPs.getOutOfRangeInd(ActrRC,IszRC);
        ActrRC(indNaN,:)=1;
        % -------------------------------------------------------------------------------


        % XYZ AT CENTER OF LE SAMPLED PATCHES IN 3D
        try
            ind=sub2ind(IszRCD,ActrRC(:,1),ActrRC(:,2));
        catch ME
            disp(['Size Axyz ' Num.toStrSane(size(Axyz)) ]);
            rethrow(ME);
        end
        ind=round([ind ind+n*m ind+2*n*m ]);

        %POINTS OUT OF RANGE
        invInd=any(ind>res,2);
        ind(invInd,:)=1; %temorparily assign 1, doesn't matter what value, kept track of by invInd

        %% NOTE
        if any(~Num.isInt(ind))
            AxyzCtr=interp2(Axyz(:,:,:),ActrRC(:,2),ActrRC(:,1));
        else
            AxyzCtr = Axyz(ind);
        end

%% Handle NaNs in AxyzCtr
        ind        = any(isnan(AxyzCtr),2);
        indNaN     = indNaN | ind | invInd;
        AxyzCtr(indNaN,:)=1;

        %%% CHECKED BELOW

        % VERGENCE ANGLE AT SAMPLED SCENE POINT
        AVrgCtrDeg = vergenceFromRangeXYZVec(A,IPDm,AxyzCtr);

        % LE AND RE LINES OF SIGHT TO LEFT EYE POINT
        A2AvctLOS = createLine3d(AxyzEye,AxyzCtr);
        A2BvctLOS = createLine3d(BxyzEye,AxyzCtr);

        % INTERSECTION OF PROJECTION PLANE (IPP) W. LINES OF SIGHT TO POINT IN L IMAGE
        A2AxyzIPP = intersectLinePlane(A2AvctLOS,PrjPln); % must be    center of pixel location
        A2BxyzIPP = intersectLinePlane(A2BvctLOS,PrjPln); % may not be center of pixel location


        % -------------------------------------------------------------------------------
        % 3d SAMPLED POINT IN IMAGE B NEAREST THE CORRESPONDING POINT (INDEX,ROW,COL)
        BctrRC=zeros(N,2);

        %ASTEST
        closestY = interp1(AppYm(:,1),AppYm(:,1),A2BxyzIPP(:,2),'nearest','extrap');
        closestX = interp1(AppXm(1,:),AppXm(1,:),A2BxyzIPP(:,1),'nearest','extrap');

        BctrRC(:,1)=ceil((closestY./distYperPix)+IszRC(1)/2);
        BctrRC(:,2)=ceil((closestX-CxyzEye(1))./distXperPix+IszRC(2)/2);

        ind=(BctrRC(:,1)==0 | BctrRC(:,2)==0);
        BctrRC(ind,1)=1;
        BctrRC(ind,2)=1;

        %XXX
        %FIND CPs OUTSIDE OF RANGE
        ind=(ActrRC(:,1) ~= BctrRC(:,1));
        if any(ind)
            %disp(['LRSIcorrespondingPointL2R: WARNING! forcing vertical pixel row to be the same.']);
            BctrRC(ind,:)=[ActrRC(ind,1) BctrRC(ind,2)];
        end
        %badInd=find(BctrRC(:,2)>IszRC(2) | BctrRC(:,2)<=0);
        indNaN=indNaN | CPs.getOutOfRangeInd(BctrRC,IszRC);

%% NOTE
        % CORRESPONDING POINTS MUST HAVE SAME VERTICAL VALUE (i.e. THEY LIE IN AN EPIPOLAR PLANE)
        %d=abs(ActrRC(:,1)-BctrRC(:,1));
        %badCol=( d <= 15 & d > 0 );
        %badProp=sum(badCol)/size(ActrRC,1);
        %if badProp > .1

        %    Error.warnSoft(['LRSIcorrespondingPointA2B: Significant number of columns being foreced to be the same. ' num2str(badProp)]);
        %end
        %BctrRC(badCol,:)=[ActrRC(badCol,1) BctrRC(badCol,2)];
        %badCol=d > 101;
        %if sum(badCol)/size(ActrRC,1) > .25
        %    Error.warnSoft('Computing CPS: Something may be wrong. Too many vertical locations misaligned');
        %    disp(sum(badCol)/size(ActrRC,1));
        %end
        %
        %indBadVert=CPs.getOutOfRangeInd(BctrRC,IszRC) | badCol;
        %BctrRC(indBadVert,:)=1;
        %indNaN=indNaN | indBadVert;
%%

        % XYZ OF RE SAMPLED POINT NEAREST THE TRUE CORRESPONDING POINT IN PROJECTION PLANE (in left coordinate frame)
        M=bsxfun(@plus,Bxyz,reshape(BxyzEye,1,1,3));
        BxyzCtr=zeros(N,3);
        sub=[repelem(BctrRC(:,1),3,1),...
             repelem(BctrRC(:,2),3,1),...
             repmat([1;2;3] ,size(BctrRC,1),1)...
            ];
        try
            ind=sub2ind(size(M),sub(:,1),sub(:,2),sub(:,3));
        catch ME
            disp(size(M));
            disp(Num.minMax(BctrRC(:,1)));
            disp(Num.minMax(BctrRC(:,2)));
            disp(rethrow(ME));
        end

        BxyzCtr = M(ind);
        BxyzCtr = transpose(reshape(BxyzCtr,3,size(BxyzCtr,1)/3));


        % VERGENCE ANGLE AT SAMPLED SCENE POINT FOR OPPOSITE EYE
        BVrgCtrDeg   = vergenceFromRangeXYZVec(A,IPDm,reshape(BxyzCtr,[N 1 3]));

        % LE AND RE LINES OF SIGHT TO RIGHT EYE POINT
        % R2LvctLOS    = createLine3d(LxyzEye,BxyzCtr);
        % R2RvctLOS    = createLine3d(RxyzEye,BxyzCtr);

        % INTERSECTION OF PROJECTION PLANE (IPP) W. LINES OF SIGHT TO POINT IN R IMAGE
        %B2AxyzIPP    = intersectLinePlane(createLine3d(AxyzEye,BxyzCtr),PrjPln); % must be    center of pixel location
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
end
end
