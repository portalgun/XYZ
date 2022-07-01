classdef XYZ_project < handle
methods
    function obj=get_itpXYZ_all_bi(obj)
        obj.get_CPs_all_bi(); % XXX if empty
        obj.get_itpXYZ_bi();
        % XXX
        %obj.imagesc_itpXYZ();
    end
    function obj=get_itpXYZ_bi(obj)
        obj.itpXYZ{1}=obj.get_itpXYZ(obj.CPs{1}{1},obj.CPs{1}{2});
        obj.itpXYZ{2}=obj.get_itpXYZ(obj.CPs{2}{1},obj.CPs{2}{2});
    end
    function itpXYZ=get_itpXYZ(obj,LitpRC,RitpRC)
        CppXm=obj.db.IppXm{3};
        CppYm=obj.db.IppYm{3};

        % CONVERT CORRESPONDING POINTS TO METERS IN CYCLOPIAN XY
        LitpXYZ = [interp1(CppXm(1,:),LitpRC(:,2)) interp1(transpose(CppYm(:,1)), LitpRC(:,1)) obj.db.IppZm.*ones(size(LitpRC,1),1)];
        RitpXYZ = [interp1(CppXm(1,:),RitpRC(:,2)) interp1(transpose(CppYm(:,1)), RitpRC(:,1)) obj.db.IppZm.*ones(size(LitpRC,1),1)];

        % FIND CYCLOPIAN SCENE COORDINATES OF POINTS BY TRACING FROM EYES THROUGH CORRESPONDING POINTS
        itpXYZ = intersectLinesFromPoints(obj.db.LExyz,LitpXYZ,obj.db.RExyz,RitpXYZ);
        itpXYZ = reshape(itpXYZ,obj.db.IszRC(2),obj.db.IszRC(1),3);
        itpXYZ = permute(itpXYZ,[2 1 3]);
    end
end
methods(Static=true)
    function pointsXYZ=forward_project(LExyz,RExyz,LitpXY,RitpXY,CppXm,CppYm,CppZm,X,Y,DppZm)
        DimgSepXm=0; % Display plane separation
        K=1; % image scaling

        % pixels to scen points
        %n=size(PPxyL,1);
        %LExyz=repmat(LExyz,n,1);
        %RExyz=repmat(RExyz,n,1);
        %Z=repmat(IppZm,n,1);
        LitpRC=fliplr(LitpXY);
        RitpRC=fliplr(RitpXY);

        if ~exist('DppZm','var') || isepmty(DppZm)
            DppZm=CppZm;
        end

        %LitpXYm = [interp1(CppXm(1,:),LitpRC(:,2)) interp1(CppYm(:,1)',LitpRC(:,1))];
        %RitpXYm = [interp1(CppXm(1,:),RitpRC(:,2)) interp1(CppYm(:,1)',RitpRC(:,1))];

        x=(1:length(CppXm(1,:)))';
        y=(1:length(CppYm(:,1)))';
        LitpXYm = [XYZ_project.interp1qr(x,CppXm(1,:),LitpRC(:,2)) XYZ_project.interp1qr(y,CppYm(:,1)',LitpRC(:,1))];
        RitpXYm = [XYZ_project.interp1qr(x,CppXm(1,:),RitpRC(:,2)) XYZ_project.interp1qr(y,CppYm(:,1)',RitpRC(:,1))];

        %Image CPs FROM OBSERVER'S PERSPECTIVE
        LitpCxyz=[LitpXYm(:,1), LitpXYm(:,2), CppZm*ones(size(LitpXYm,1),1)];
        RitpCxyz=[RitpXYm(:,1), RitpXYm(:,2), CppZm*ones(size(RitpXYm,1),1)];

        % meters to xyz
        LitpDxyz=[K*LitpCxyz(:,1)-DimgSepXm/2, K*LitpCxyz(:,2), DppZm*ones(size(LitpCxyz,1),1)];
        RitpDxyz=[K*RitpCxyz(:,1)+DimgSepXm/2, K*RitpCxyz(:,2), DppZm*ones(size(LitpCxyz,1),1)];

        pointsXYZ=intersectLinesFromPoints(LExyz,LitpDxyz,RExyz,RitpDxyz);
    end
    function [LitpRC,RitpRC]=back_project(LExyz,RExyz,pointsXYZ,IppXm,IppYm,IppZm,X,Y)
    %% scene points to pixels
        n=size(pointsXYZ,1);
        LExyz=repmat(LExyz,n,1);
        RExyz=repmat(RExyz,n,1);

        Lline=createLine3d(LExyz,pointsXYZ);
        Rline=createLine3d(RExyz,pointsXYZ);
        plane=createPlane([1 0 IppZm],[0 1 IppZm],[0 -1 IppZm]);

        PPxyzLM=intersectLinePlane(Lline,plane);
        PPxyzRM=intersectLinePlane(Rline,plane);

        LitpRC=zeros(size(PPxyzLM,1),2);
        RitpRC=zeros(size(PPxyzRM,1),2);

        [LitpRC(:,2),LitpRC(:,1)]=XYZ.revImgInterp2(IppXm,IppYm,X,Y,PPxyzLM(:,1),PPxyzLM(:,2));
        [RitpRC(:,2),RitpRC(:,1)]=XYZ.revImgInterp2(IppXm,IppYm,X,Y,PPxyzRM(:,1),PPxyzRM(:,2));

    end
    function  [pnL,pnR]=interpScattered(PPxyL,PPxyR,IszRC,map)
        pmap=(map(:));
        [X,Y]=meshgrid(1:IszRC(2),1:IszRC(1));

        indNanL=any(isnan(PPxyL),2);
        indNanR=any(isnan(PPxyR),2);
        F=scatteredInterpolant(PPxyL(~indNanL,1),PPxyL(~indNanL,2),pmap(~indNanL),'natural','none');
        pnL=F(X,Y);


        F=scatteredInterpolant(PPxyR(~indNanR,1),PPxyR(~indNanR,2),pmap(~indNanR),'natural','none');
        pnR=F(X,Y);
    end
    function xyzM=interp(PPxy,IppXm,IppYm)
        xyzM=Map.interp2(IppXm,IppYm,X,Y,PPxy(:,1),PPxy(:,2));
    end
    %function PPxy=revInterp(xyzM,IppXm,IppYm,X,Y)
    %    PPxy=Map.revImgInterp2(IppXm,IppYm,X,Y,xyzM(:,1),xyzM(:,2));
    %end
    function [Xitp,Yitp]=revImgInterp2(Xm,Ym,X,Y,Vxm,Vym)
        % XXX CHECK MINUS -0.5
        %ind=size(X,2)/2
        %X(1,ind)
        %X(1,ind+1)
        add=-0.5;
        if Arr.ndim(Vxm)==1
            %Xitp=interp1(Xm(1,:),X(1,:)+add,Vxm,'linear');
            %Yitp1=interp1(Ym(:,1),Y(:,1)+add,Vym,'linear');

            Xitp=XYZ_project.interp1qr(Xm(1,:)',X(1,:)'+add,Vxm);
            Yitp=XYZ_project.interp1qr(flip(Ym(:,1),1),flip(Y(:,1),1)+add,flip(Vym,1));
        else
            Xitp=interp2(Xm,Ym,X(:)+add,Vxm,Vym,'linear');
            Yitp=interp2(Xm,Ym,Y(:)+add,Vxm,Vym,'linear');
        end
    end
    function [Xitp,Yitp]=imgInterp2(X,Y,Xm,Ym,Vx,Vy)

    % X Y in pixels
    % Ym Xm is meters projection plane
    %
        if Arr.ndim(Vx)==1
            Xitp=interp1(X(1,:),Xm(1,:),Vx,'linear');
            Yitp=interp1(Y(:,1),Ym(:,1),Vy,'linear');
        else
            Xitp=interp2(X,Y,Xm(:),Vx,Vy,'linear');
            Yitp=interp2(X,Y,Ym(:),Vx,Vy,'linear');
        end
    end
    function yi=interp1F(x,y,xi,ydiff)
        % x is pix
        % y is m

        m = size(x,1);

        % ind1 = hist(ci,x).bin = round(xi) = x(ind1)
        ind1=floor(xi);  % XXX bottlneck 2

        ind1 = max(ind1,1);     % To avoid index=0 when xi < x(1)
        ind1 = min(ind1,m-1);   % To avoid index=m+1 when xi > x(end).

        yi = y(ind1) + (xi-ind1)*ydiff; % XXX bottlneck 1
    end
    function yi=interp1qr(x,y,xi)
        % XXX NOT USED?
    %https://www.mathworks.com/matlabcentral/fileexchange/43325-quicker-1d-linear-interpolation-interp1qr
        %x=Vec.col(x);
        m = size(x,1);
        %n = size(y,2);

        % For each 'xi', get the position of the 'x' element bounding it on the left [p x 1]
        [~,ind1] = histc(xi,x); % XXX bottleneck 1, get left bound

        %max(ind1)
        %min(ind1)
        %max(ind1)
        %m-1
        ind1 = max(ind1,1);     % To avoid index=0 when xi < x(1)
        ind1 = min(ind1,m-1);   % To avoid index=m+1 when xi > x(end).

        ind2=ind1+1;
        x1=x(ind1);
        y1=y(ind1);


        t = (xi-x1)./(x(ind2)-x1); % dxi/dx
        % Get 'yi'
        yi = y1 + t.*(y(ind2)-y1); % XXX bottlneck 2
        % Give NaN to the values of 'yi' corresponding to 'xi' out of the range of 'x'


        %yi(xi<x(1) | xi>x(end),:) = NaN;
    end


end
end
