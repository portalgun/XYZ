classdef XYZ_vrg < handle
methods(Static=true)
    function Ivrg_deg=toVrg(LorRorC,IPDm,Ixyz,dim)
        [C,LExyz,RExyz,dim,bVec]=XYZ.parse_to_vrg__(LorRorC,IPDm,Ixyz,dim);
        if bVec
            Ivrg_deg=XYZ.to_vrg_vec_(C.^2,LExyz,RExyz,Ixyz,dim);
        else
            Ivrg_deg=XYZ.to_vrg_(C,LExyz,RExyz,Ixyz,dim);
        end
    end
    function IvrgDeg=list_to_vrg_vec_(Csqr,LExyz,RExyz,Ixyz,dim)
        if isequal(dim,0)
            IvrgDeg=XYZ.to_vrg_naive(RExyz(1)-LExyz(1),Ixyz(:,3));
            return
        end
        Asqr=sum( (Ixyz(:,[dim 3])-RExyz).^2,2);
        Bsqr=sum( (Ixyz(:,[dim 3])-LExyz).^2,2);

        IvrgDeg = acosd( (Asqr + Bsqr - Csqr)./(2.*realsqrt(Asqr.*Bsqr)));
        negInd=Ixyz(:,3)<=0;
        IvrgDeg(negInd)=IvrgDeg(negInd)*-1;
    end
    function IvrgDeg=to_vrg_vec_(Csqr,LExyz,RExyz,Ixyz,dim)
        if isequal(dim,0)
            IvrgDeg=XYZ.to_vrg_naive(RExyz(1)-LExyz(1),Ixyz(:,:,3));
            return
        end
        Asqr=sum( (Ixyz(:,:,[dim 3])-RExyz).^2,3);
        Bsqr=sum( (Ixyz(:,:,[dim 3])-LExyz).^2,3);

        IvrgDeg = acosd( (Asqr + Bsqr - Csqr)./(2.*realsqrt(Asqr.*Bsqr)));
        negInd=Ixyz(:,:,3)<=0;
        IvrgDeg(negInd)=IvrgDeg(negInd)*-1;
    end

    function IvrgDeg=to_vrg_(C,LExyz,RExyz,Ixyz,dim)
        if isequal(dim,0)
            IvrgDeg=XYZ.to_vrg_naive(RExyz(1)-LExyz(1),Ixyz(:,:,3));
            return
        end

        Asqr=sum(bsxfun(@minus,Ixyz(:,:,[dim 3]),RExyz).^2,3);
        Bsqr=sum(bsxfun(@minus,Ixyz(:,:,[dim 3]),LExyz).^2,3);

        IvrgDeg = acosd( (Asqr + Bsqr - C.^2)./(2.*realsqrt(Asqr.*Bsqr)));

        %IvrgDeg = XYZ.SSSdeg__(A,B,C);
        negInd=Ixyz(:,:,3)<=0;
        IvrgDeg(negInd)=IvrgDeg(negInd)*-1;
    end
    function IvrgDeg=to_vrg_naive(IPDm,z);
        % Z coordinate in meters of query point
        IvrgDeg=2.*atand(IPDm./(2.*z));
    end
    function thetaDeg=SSSdeg__(a,b,c)
        thetaDeg = acosd( (a.^2 + b.^2 - c.^2)./(2.*a.*b));
    end
    function [C,LExyz,RExyz,dim,bVec] = parse_to_vrg__(LorRorC,IPDm,IxyzORIszRCD,dim,bVec)
        % TODO IPDm to LExyz Rxyz

        if ~exist('bVec','var') || isempty(bVec)
            bVec=false;
        end

        if numel(IxyzORIszRCD)==3
            IszRCD=IxyzORIszRCD;
        else
            Ixyz=IxyzORIszRCD;
            if size(Ixyz,3) ~= 3
                error(['vergenceFromRangeXYZ: WARNING! size(Ixyz,3) = ' num2str(size(Ixyz,4)) ' instead of 3']);
            end
        end
        if ~exist('dim','var') || isempty(dim)
            dim=[1:2];
        end

        switch LorRorC
        case 'L'
            LExyz = [   0 0 0];
            RExyz = [IPDm 0 0];
        case 'R'
            LExyz = [-IPDm 0 0];
            RExyz = [    0 0 0];

        case 'C'
            LExyz = [-IPDm/2 0 0];
            RExyz = [ IPDm/2 0 0];
        otherwise
            error(['vergenceAngleFromRangeXYZ: WARNING! unhandled LorRorC string: ' LorRorC]);
        end
        if numel(dim)==1 & dim==2
            C=0;
        else
            C=IPDm;
        end
        if numel(dim)==1 & dim==0
           ;
        else
            sz=[1,1,numel(dim)+1];
            LExyz=reshape(LExyz([dim 3]),sz);
            RExyz=reshape(RExyz([dim 3]),sz);
            if bVec
                r=size(Ixyz(:,:,1),1);
                c=size(Ixyz(:,:,1),2);
                LExyz=repmat(LExyz,r,c,1);
                RExyz=repmat(RExyz,r,c,1);
            end
        end

    end
    function [vrg,vrsC,vrsL,vrsR]=toVrgAngle(xyz,LExyz,RExyz)
        % TODO - check with toVrg
        if nargin < 2
            LExyz=[-0.065/2 0 0];
            RExyz=[ 0.065/2 0 0];
        end

        LExy=LExyz(1:2);
        RExy=RExyz(1:2);
        bReshape=false;
        if numel(xyz)==3
            xy=xyz(1:2);
            z=xyz(end);
        elseif size(xyz,3)==1 && size(xyz,2) == 3
            xy=xyz(:,1:2);
            z=xyz(:,3);
        elseif size(xyz,1)==size(xyz,2) && size(xyz,3)==3
            x=xyz(:,:,1);
            y=xyz(:,:,2);
            z=xyz(:,:,3);
            xy=[x(:) y(:)];
            z=z(:);
        else
            error('xyz not of proper dimensions')
        end
        Ez=LExyz(3);
        z=z-Ez;

        vrsL=atand((xy-LExy)./z);
        vrsR=atand((xy-RExy)./z);
        vrg=vrsL - vrsR;
        if nargout > 1
            vrsC=atand(xy./z);
        end

        %out(1)=vrg_fun(xyz(1),xyz(3),LExyz(1),RExyz(1));
        %out(2)=vrg_fun(xyz(2),xyz(3),LExyz(2),RExyz(2));

        %function out=vrg_fun(x,z,LEx,REx)

            %if     x < 0
            %    A=abs(LEx);
            %    B=abs(REx);
            %elseif x > 0
            %    A=abs(REx);
            %    B=abs(LEx);
            %else
            %    I=REx-LEx;
            %    I
            %    out=2.*atand(I./(2.*z));
            %    return
            %end
            %out=atand((x+B)/z) - atand((x-A)/z);
        %end
    end
    function toVrgAngle_test_()
        n=100;

        x=linspace(-1,1,n);
        y=zeros(n,1);
        z=ones(n,1);
        xyz=[x' y z];

        v=XYZ.toVrgAngle(xyz);

        %v=zeros(n,2);
        %for i = 1:n
        %    v(i,:)=XYZ.get_vrg_vrs([x(i) 0 1]);
        %end

        subplot(1,2,1);
        hold off;
        plot(x,v(:,1));
        xlabel('x');
        ylabel('vrg');
        axis square;


        x=zeros(n,1);
        y=zeros(n,1);
        z=linspace(0,1,n);
        xyz=[x y z'];

        v=XYZ.toVrgAngle(xyz);

        %v=zeros(n,2);
        %for i = 1:n
        %    v(i,:)=XYZ.get_vrg_vrs([0 0 z(i)]);
        %end
        subplot(1,2,2);
        hold off;
        plot(z,v(:,1));
        xlabel('z');
        ylabel('vrg');
        axis square;
    end
    function [vrgXY,vrsXY]=get_vrg_vrs_map(xyz,LExyz,RExyz)
        % XXX BROKEN
        % XXX should I use xyz2vrgVrs code?
       % get_vrg_vrs_map(xyz, LExyz, RExyz)
       %if Arr.ndim(xyz)==2 && size(xyz,2)==3

        error('broken')


        CExyz=(LExyz+RExyz)./2;
        IPDm=sqrt(sum(RExyz-LExyz,2).^2);
        I=IPDm/2;
        %%%%
        L=sqrt(sum(xyz-LExyz,2).^2);
        R=sqrt(sum(xyz-RExyz,2).^2);
        C=sqrt(sum(xyz-CExyz,2).^2);

        vrgDeg=acosd((L.^2+R.^2-IPDm)./(2.*L.*R));

        rdeg=acosd((I.^2+R.^2-C.^2)./(2.*I.*R));
        ldeg=180-vrgDeg-rdeg;
        cdeg=asind(R.*sind(rdeg)./C);

        vrsDeg=cell(1,3);
        vrsDeg{1}=90-ldeg;
        v1=vrsDeg{1}(ldeg > 90);
        if ~isempty(v1)
            vrsDeg{1}(ldeg > 90)=v1.*-1;
        end

        vrsDeg{3}=90-cdeg;
        v3=vrsDeg{1}(cdeg > 90);
        if ~isempty(v3)
            vrsDeg{3}(cdeg > 90)=v3.*-1;
        end

        vrsDeg{2}=90-rdeg;
        v2=vrsDeg{2}(rdeg <= 90);
        if ~isempty(v2)
            vrsDeg{2}(ldeg <= 90)=v2.*-1;
        end
        %%%
    end

    function z=vrg_vrs_map_to_xyz(vrgDeg,vrsDeg,LExyz,RExyz)
        % XXX is version defined correctly here?
        IPDm=sqrt(sum(RExyz-LExyz,2).^2);
        error('unverified')

        C=acotd(2*tand(vrsDeg));
        L=-1*vrgDeg + C; %left angle
        R=180-vrgDeg-L;                       %right angle
        s=sind(L)*IPDm/sind(vrgDeg);          %solve for one side
        z=sind(R)*s;

    end
    function vrg=vrsAndDistToVrg(vrsDeg,dist)
        z=sind(90-abs(vrsDeg)).*dist;
        vrg=acosd(z/dist);
    end
    function [z,z1,z2]=vrgAndXToZ(vrgDeg,xy,IPDm,dim)
        % inverse of toVrg
        % convert vergence & x coordinate to z coordinate
        sz=size(vrgDeg);
        I=IPDm;
        g=(vrgDeg(:)*pi/180);

        if dim==1
            X=xy(:,:,1);
        elseif dim == 2
            X=xy(:,:,2);
        else
            error(['dimensions must be 1 or 2'])
        end
        X=(X(:));

        [z1,z2]=inv_fun(I,g,X);

        z=nan(length(g),1);
        ind1=g>=0 & z1>=0 & z2<0;
        z(ind1)=z1(ind1);

        ind1=g<0 & z1<0 & z2>=0;
        z(ind1)=z1(ind1);


        ind2=g>=0 & z2>=0 & z1<0;
        z(ind2)=z1(ind2);

        ind1=g<0 & z2<0 & z1>=0;
        z(ind2)=z1(ind2);

        ind=z1==z2;
        z(ind)=z1(ind);

        z=reshape(z,sz);
        z1=reshape(z1,sz);
        z2=reshape(z2,sz);

        z(isnan(z))=0;

        function [z1,z2]=inv_fun(I,g,x)
            nInd=ones(size(x));
            nInd((abs(x) <= IPDm/2))=-1;

            h(:,1)=nInd.* atan((I - (I.^2.*tan(g).^2 - 4.*x.^2.*tan(g).^2 + I.^2).^(1/2))./(I.*tan(g) + 2*x.*tan(g)));
            h(:,2)=nInd.* atan((I + (I.^2.*tan(g).^2 - 4.*x.^2.*tan(g).^2 + I.^2).^(1/2))./(I.*tan(g) + 2*x.*tan(g)));

            z2=(x+I/2)./tan(g+nInd.*h(:,1));
            z1=(x+I/2)./tan(g+nInd.*h(:,2));
        end

    end
    function [vrsL,vrsR]=get_eye_versions(obj,LExyz,RExyz,posXYZm)
        % XXX check
        straightL=[Lxyz(1) LExyz(2) posXYZm(3)];
        straightR=[Rxyz(1) RExyz(2) posXYZm(3)];
        L=xyz2triangleAngles(LExyz,straightL,posXYZm);
        R=xyz2triangleAngles(RExyz,straightR,posXYZm);
        vrsL=rad2deg(L);
        vrsR=rad2deg(R);
    end
    function elev=get_elev(obj,trgtXYZ,posXYZm)
        % XXX
        elev=atand(trgtXYZ(2)/posXYZm(3));      %ELEVATION ...
    end
end
end
