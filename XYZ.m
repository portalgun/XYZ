classdef XYZ < handle & dbImg & idxConversions & XYZ_plot & XYZ_project & XYZ_vrg & XYZ_disparity
%xyz=XYZ('LRSI',1)
%xyz.get_CPs_all_bi
% LorRorC:      left,right or cyclopean co-ordinate system to use
%            'L' -> left      image used   as reference
%            'R' -> right     image used   as reference
%          ? 'C' -> cyclopean image ?
% IPDm:      inter-ocular distance to simulate
% Ixyz:      range data in cartesian coordinates [ r x c x 3 ]
%
%            Ixyz(:,:,1) -> x-values
%            Ixyz(:,:,2) -> y-values
%            Ixyz(:,:,3) -> z-values
properties
    xyz
    itpXYZ
    dispXYZ

    CPs
    LitpRC
    RitpRC
    CitpRC
    LitpRCchk
    RitpRCchk

    % db
    disp

    IvrgDeg
end
properties(Hidden=true)
    MLR
    cpLookup
end
methods
    function obj=XYZ(database,I, db,MLR,allRC,allXYZ)
        if ~exist('db','var')
            db=[];
        end
        obj@dbImg(database,'img','xyz',I,[],0,db);
        obj.xyz{1}=obj.im.xyz{1};
        obj.xyz{2}=obj.im.xyz{2};
        obj.im=[];

        if exist('MLR','var') && ~isempty(MLR)
            obj.MLR=MLR;
        end

        if exist('allRC','var') && ~isempty(allRC)
            obj.allRC=allRC;
        end
        if exist('allXYZ','var') && ~isempty(allXYZ)
            obj.allXYZ=allXYZ;
        end

    end
    function obj=set_image(obj,I)
        % XXX ?
        obj.I=I;
        obj.reset_CPs();
        obj.reset_cpLookup();
    end
end
methods(Static=true)
%%- FILES
%% BASE
    function name=getName(I,k)
        if ischar(i)
            LorR=k;
        elseif Num.is(k)
            LorR=CPs.LANDR(k);;
        end
        name=[LorR num2str(I,'%03i')];
    end
    function fname=getFname(database,I,K)
        dire=XYZ.getDire(database);
        name=XYZ.getName(I,K);
        fname=[dire name];
    end
    function dire=getDire(database)
        dire=Env.var('ImgDb.xyz',database);
    end
    function xyz=load(database,I,k)
        fname=XYZ.getFName(database,I,k);
        load(fname);
    end
%% NaN
    function name=getNameNan(I,k)
        if ischar(i)
            LorR=k;
        elseif Num.is(k)
            LorR=CPs.LANDR(k);;
        end
        name=[LorR num2str(I,'%03i')];
    end
    function fname=getFnameNan(database,I,K)
        dire=XYZ.getDireNan(database);
        name=XYZ.getNameNan(I,K);
        fname=[dire name];
    end
    function dire=getDireNan(database)
        dire=Env.var('ImgDb.xyzNan',database);
    end
    function xyz=loadNan(database,I,k)
        fname=XYZ.getFNameNan(database,I,k);
        load(fname);
    end
%% InPaint
    function name=getNameInPaint(I,k)
        if ischar(i)
            LorR=k;
        elseif Num.is(k)
            LorR=CPs.LANDR(k);;
        end
        name=[LorR num2str(I,'%03i')];
    end
    function fname=getFnameInPaint(database,I,K)
        dire=XYZ.getDireInPaint(database);
        name=XYZ.getNameInPaint(I,K);
        fname=[dire name];
    end
    function dire=getDireInPaint(database)
        dire=Env.var('ImgDb.xyzInPaint',database);
    end
    function xyz=loadInPaint(database,I,k)
        fname=XYZ.getFNameInPaint(database,I,k);
        load(fname);
    end
%%- Misc
    function z=elevationToZ(elev,y)
        z=y./tand(elev);
    end
    function out=getWHD(xyz)
        out=diff( [min(xyz,[],[1,2]) max(xyz,[],[1,2])] ); % SLOW
        out=abs(out(:));
    end
%%- ANGLES
    function theta=angle(p1xyz,p2xyz,p3xyz)
        A=p1xyz;
        B=p2xyz;
        C=p3xyz;
        v1=A-B;
        v2=C-B;


        v1mag=sqrt(sum(v1.^2));
        v2mag=sqrt(sum(v2.^2));
        v1norm=v1/v1mag;
        v2norm=v2/v2mag;
        res=dot(v1norm,v2norm');
        theta = acos(res);

    end
    function theta=angles(M,N,O)
        for i = 1:3
            j=i+1;
            if j == 4
                j=1;
            end
            m=M; n=N; o=O;
            m(i)=0; n(i)=0; o(i)=0;
            theta(i)=XYZ.angle(m,n,o);
        end
        theta=[-theta(1) theta(2) theta(3)];
        %theta=[theta(1) theta(3) theta(2)];
        %theta=[theta(2) theta(1) theta(3)];
        %theta=[theta(2) theta(3) theta(1)];
        %theta=[theta(3) theta(1) theta(2)];
        %theta=[theta(3) theta(2) theta(1)];
    end
    function elev=elevation(xyz)
        % 0-90 = above & ahead
        % 90-180 = above & behind
        % -0--90 = below & ahead
        % -90--180 = below & behind
        y=xyz(:,:,2);
        z=xyz(:,:,3);
        elev=atand(y./z);
        below=y<=0;
        behind=z<=0;
        above=~below;
        front=~behind;

        % + above & front   1 0-90       0-90
        % - above & behind  2 0-90     90-180
        % + below & behind  3 -90-90   -270--90
        % - below & front   4 -90-0    =-0--90

        %elev(above & front) =nan;
        elev(above & front) =elev(above & front);
        %Num.minMax(elev(above & front))
        %
        %elev(above & behind)=nan;
        elev(above & behind)=elev(above & behind)+180;
        %Num.minMax(elev(above & behind))

        %elev(behind & below)=nan
        elev(below & behind)=elev(below & behind)-180;
        %elev(below & behind)=elev(below & behind)-0;
        %Num.minMax(elev(below & behind))

        %elev(below & front) =nan
        elev(below & front) =elev(below & front);
        %Num.minMax(elev(below & front))

    end
%%- WARP
    function xyz=transform_to_win(xyz,winPosXYZm,winWHm,PszRC)
        % XXX RENAME
        xyz=XYZ.scale_and_center(xyz,winWHm,PszRC,winPosXYZm);
        xyz=XYZ.flatSheer(xyz,winWHm);   %% bottleneck 2
    end
    function xyz=scale_and_center(xyz,WHm,PszRC,posXYZm)
        % XXX RENAME
        PszRCCur=size(xyz);
        ctr=floor(PszRCCur/2+1); % verified

        %mult(3)=1;

        xyz=XYZ.scale_and_center_(xyz,ctr,mult,posXYZm);
    end
    function xyz=rotate(obj,xyz,xyzTheta)
        % TEST
        x=xyz(:,:,1);
        y=xyz(:,:,2);
        z=xyz(:,:,3);

        xyz(:,:,2) = y*cos(xyzTheta) - z*sin(xyzTheta);
        xyz(:,:,3) = y*sin(xyzTheta) + z*cos(xyzTheta);
        y=xyz(:,:,2);
        z=xyz(:,:,3);

        xyz(:,:,1) = x*cos(xyzTheta) + z*sin(xyzTheta);
        xyz(:,:,3) = z*cos(xyzTheta) - x*sin(xyzTheta);
        x=xyz(:,:,1);
        z=xyz(:,:,3);

        xyz(:,:,1) = x*cos(xyzTheta) - y*sin(xyzTheta);
        xyz(:,:,2) = x*sin(xyzTheta) + y*cos(xyzTheta);
    end
    function xyz=scale(xyz,WHm,PszRC)
        PszRCCur=size(xyz);
        ctr=floor(PszRCCur/2+1); % verified

        buffmult=PszRCCur(1:2)./PszRC;

        WHD=XYZ.getWHD(xyz);
        ZtoY=WHD(3)/WHD(2);
        mult(1)=WHD(1)/(WHm(1)*buffmult(2));
        mult(2)=WHD(2)/(WHm(2)*buffmult(1));
        mult(3)=ZtoY*WHD(2)/(WHm(2));
        %mult(3)=1;

        xyz=XYZ.scale_(xyz,ctr,mult);


    end
    function xyz=center(xyz,posXYZm)
        ctr=floor(size(xyz)/2+1); % verified
        for i = 1:3
            % center at zero & center at posXYZM
            xyz(:,:,i)=xyz(:,:,i)-xyz(ctr(1),ctr(2),i) + posXYZm(i);;
        end
    end
    function xyz=flatSheer(xyz,WHm);
       PszRC=size(xyz);
       W=WHm(1)/2;
       H=WHm(2)/2;
       [X,Y]=meshgrid(linspace(-W,W,PszRC(2)),linspace(H,-H,PszRC(1)));
       xyz(:,:,1)=X;
       xyz(:,:,2)=Y;
    end
end
methods(Static, Access=?XYZ_disparity)
    function mult=get_scale_mult__(xyz,PszRC,PszRCCur);
        buffmult=PszRCCur(1:2)./PszRC;

        WHD=XYZ.getWHD(xyz);
        ZtoY=WHD(3)/WHD(2);
        mult(1)=WHD(1)/(WHm(1)*buffmult(2));
        mult(2)=WHD(2)/(WHm(2)*buffmult(1));
        mult(3)=ZtoY*WHD(2)/(WHm(2));
    end
    function mult=get_scale_mult_z__(xyz,WHm,PszRC,PszRCCur);
        buffmult=PszRCCur(1:2)./PszRC;

        WHD=XYZ.getWHD(xyz);
        ZtoY=WHD(3)/WHD(2);
        mult=ZtoY*WHD(2)/(WHm(2));
    end
end
methods(Static, Hidden)
    function xyz=flatSheer_(xyz,Xm,Ym);
       xyz(:,:,1)=Xm;
       xyz(:,:,2)=Ym;
    end
    function xyz=scale_(xyz,ctr,mult)


        for i = 1:3
            c=xyz(ctr(1),ctr(2),i);

            % center @ 0, scale, recenter
            xyz(:,:,i)=((xyz(:,:,i)-c)./mult(i)) + c;
        end
    end
    function xyz=scale_and_center_(xyz,ctr,mult,posXYZm)
        for i = 1:3
            c=xyz(ctr(1),ctr(2),i);
            % center at zero & center at posXYZM
            xyz(:,:,i)=((xyz(:,:,i)-c)./mult(i)) + posXYZm(i); % SLOW
        end
    end
    function xyz=scale_center_sheer_(xyz,ctr,mult,posXYZm,X,Y)
        c=xyz(ctr(1),ctr(2),3);
        xyz(:,:,1)=X;
        xyz(:,:,2)=Y;
        xyz(:,:,3)=((xyz(:,:,3)-c)./mult) + posXYZm(3); % SLOW
    end
end
end
