classdef XYZ_disparity < handle
methods(Static)
    function im=dnSmp(I,k,varargin)
        [N,M]=size(I);
        n=(N/k);
        m=(M/k);
        %k2=(k*[n+k m+k])./[n m];
        n2=ceil((N+k)/k);
        m2=ceil((M+k)/k);

        %122.333
        %163.333
        im=cell(k,k);
        for j = -1:k-2
        for i = -1:k-2
            %ii=k-i;
            %jj=k-j;
            trans=imtranslate(I,[-j,-i]);
            im{i+2,j+2}=imresize(trans,[n,m],varargin{:});

            %size(Itmp)
            %Itmp=[[r; Itmp] c];
            %size(Itmp)
            %Itmp=imresize(Itmp,[n2 m2]);
            %dk
        end
        end
    end
    function IM=dnSmp_test(k)
        method='nearest';
        PszRC=[51 51];
        IMo=Msk.circle(PszRC,20,[0 0],1);
        IMo=IMo([1:25 1:26],[26:51 1:25]);
        PszDnk=PszRC/k;
        if k > 1
            Im=XYZ_disparity.dnSmp(IMo,k,'Method',method);
            [KX,KY]=meshgrid(repmat(1:k,1,PszDnk(2)),repmat(1:k,1,PszDnk(1)));
            [X,Y]  =meshgrid(repelem(1:PszDnk(2),1,k),repelem(1:PszDnk(2),1,k));
        end
        %kR=arrayfun(@(x) find(mod(x,k:-1:1,1,'first')),1:PszRC(1));
        %kC=arrayfun(@(x) find(mod(x,k:-1:1,1,'first')),1:PszRC(2));
        IM=nan(PszRC);
        I=0;
        J=1;
        for j = 1:PszRC(2)
        for i = 1:PszRC(1)
            if k > 1
                im=Im{KY(i,j),KX(i,j)};
                I=Y(i,j);
                J=X(i,j);
            else
                im=Im;
                I=i;
                J=j;
            end
            try
                IM(i,j)=im(I,J);
            catch ME
                %[I J]
                %size(im)
                %rethrow(ME)
            end
        end
        end
        N=1;
        M=3;
        figure(1)
        subPlot(N,M,1,1);
        imagesc(IMo);
        title('original');
        axis square;

        subPlot(N,M,1,2);
        imagesc(IM);
        title('reconstructed');
        axis square;

        subPlot(N,M,1,3);
        imagesc(IMo-(IM>=0.5));
        title('diff');
        axis square;
        colorbar;

        % subPlot(N,M,1,3);
        % imagesc(KX);
        % axis square;

        % subPlot(N,M,1,4);
        % imagesc(KX);
        % axis square;

        % subPlot(N,M,1,5);
        % imagesc(Y);
        % axis square;

        % subPlot(N,M,1,6);
        % imagesc(X);
        % axis square;

        figure(2)
        for i = 1:k
        for j = 1:k
            subPlot(k,k,i,j);
            imagesc(Im{i,j});
        end
        end

    end
    function [Ximg,Ximg2]=test_by_XYZ()
        % XXX interpolated nan images - xyzNan
        xyz=dbImg.getImg('LRSI','img','xyzNan',1,'L');
        kernSz=[32 32];
        winPosXYZm=[0 0 1];
        winWHdeg=[1 1];

        hostname='jburge-wheatstone@psych.upenn.edu';
        vDisp=VDisp(hostname);
        winWHm=winWHdeg./vDisp.degPerMxy;
        %winWHm=[0.0178 0.0178];

        dim=1;
        dnk=4;

        %dim=1;
        [Ximg,Ximg2]=XYZ.disparity_contrast_img(xyz,kernSz,winPosXYZm,winWHm,'dim',dim,'dnk',dnk);
    end

    function [Ximg]=depth_contrast_by_CPs(cpLookup,IszRC,kernSz,k,winWHm,winPosXYZm,vDisp,varargin)
        %DB:    LExyz,RExyz,IszRC,cpLookup,
        %VDisp: IppXm,IppZm,Xpix,Yix
        %Param: PszRC,k,

        kernSz=kernSz+double(mod(kernSz,2)==0);
        S=struct(varargin{:});
        if isfield(S,'dnk') && S.dnk ~=1
            kernSzDnk=round(kernSz/S.dnk);
            kernSzDnk=kernSzDnk+double(mod(kernSzDnk,2)==0);
        else
            kernSzDnk=kernSz;
        end
        Opts=XYZ.parse_disparity_contrast__(kernSzDnk,varargin{:});
        % ENSURE ODD kernSz

        Opts.bConv=false;

        % M MESH
        W=winWHm(1)/2;
        H=winWHm(2)/2;

        % FORWARD INTERP
        ptchNumel=prod(kernSz);
        LExyzVec=repmat(vDisp.LExyz,ptchNumel,1);
        RExyzVec=repmat(vDisp.RExyz,ptchNumel,1);
        Zvec=repmat(vDisp.Zm,ptchNumel,1);
        cinit=zeros(size(RExyzVec));

        % WINDOW
        WW=Opts.W(:);
        wSum=sum(WW);

        % ITER
        Ximg=nan(IszRC);
        ctrRC=floor(kernSz./2+1);
        ctr=sub2ind(kernSz,ctrRC(1),ctrRC(2));
        Rh=floor(kernSz(1)/2);
        Rw=floor(kernSz(2)/2);
        ROWS = Rh+1:IszRC(1)-Rh-1;
        COLS = Rw+1:IszRC(2)-Rw-1;

        % Prog
        N=numel(ROWS)*numel(COLS);
        modIter=floor(N/100);
        %pr=Pr(N,modIter,'disparity contrast img');
        % 5 min per patch -> 7 hours for set


        % WINDOW M TO PIX
        % XXX check XYZ back/forward project
        WszRCpix=fliplr(winWHm.*vDisp.pixPerMxy);
        [winPosRCpixL,winPosRCpixR] =vDisp.PP.back_project(winPosXYZm,true);

        AitpRC0=CPs.getAitpRC0(kernSz);


        % PATCH TO WIN SCALE
        %cpLookup=cpLookup{k};

        for i = ROWS
        for j = COLS
            ActrRC=[i,j];


            % scene to floating patch
            [LitpRC,RitpRC,LctrRC,RctrRC]=CPs.getPatchFast(ActrRC,AitpRC0,k,cpLookup,IszRC,true); % XXX Bottleneck 2 (12%)
            if isempty(LitpRC) || any(isnan(LitpRC),'all') || any(isnan(RitpRC),'all')
                continue
            end

            % float patch to centered disp.win
            [LitpRCw,RitpRCw]=CPs.patchPixToWinPix(LitpRC,RitpRC,WszRCpix,kernSz,winPosRCpixL,winPosRCpixR);

            % centered win to XYZ
            xyz=vDisp.PP.forward_project(LitpRCw,RitpRCw, LExyzVec,RExyzVec,Zvec,cinit); % XXX Bottleneck 1 (74%) : interp, intersectLinesFromPoints


            % DOWNSAMPLE HERE
            if Opts.dnk~=1
                xyz=reshape(xyz,[kernSz 3]);
                z=imresize(xyz(:,:,3),kernSzDnk,'bilinear');


                DC   = sum(z(:).*WW)./wSum;
                Ximg(i,j)  = single(real(sqrt( sum(WW.*( z(:)-DC ).^2)./wSum )));
            else
                DC   = sum(xyz(:,3).*WW)./wSum;
                Ximg(i,j)  = single(real(sqrt( sum(WW.*( xyz(:,3)-DC ).^2)./wSum )));
            end



            bPlot=false; %NOTE
            if bPlot
                disp(['ActrRC ' Num.toStr(ActrRC) newline ...
                      'LctrRC ' Num.toStr(LctrRC) newline ...
                      'RctrRC ' Num.toStr(RctrRC) newline ]);

                Fig.newU(1);
                fr=4;

                subplot(fr,3,1)
                Plot.RC3(xyz,'.k');
                colormap hot;
                title('XYZm');
                ylim([.85 1.15]);
                axis square;

                subplot(fr,3,2)
                imagesc(reshape(vrgImg,Opts.PszRC)');
                set(gca,'YDir','normal');
                Fig.formatIm;
                colormap hot;
                colorbar;
                Fig.format('','','vrg deg');

                subplot(fr,3,3)
                imagesc(reshape(dsp,Opts.PszRC)');
                set(gca,'YDir','normal');
                colorbar;
                Fig.formatIm;
                colormap hot;
                Fig.format('','','dsp deg');

                % ---

                subplot(fr,3,4)
                hist(xyz(:,1));
                title('x');

                subplot(fr,3,5)
                hist(xyz(:,2));
                title('y');

                subplot(fr,3,6)
                hist(xyz(:,3));
                title('z');

                % ---

                subplot(fr,3,7)
                hist(LitpRC(:,1));
                title('LitpR');

                subplot(fr,3,8)
                hist(RitpRC(:,1));
                title('RitpR');

                subplot(fr,3,9)
                hist(RitpRC(:,2));
                title('RitpC');


                % ---

                subplot(fr,3,10)
                hist(LitpRCw(:,1));
                title('LitpRw');

                subplot(fr,3,11)
                hist(RitpRCw(:,1));
                title('RitpRw');

                subplot(fr,3,12)
                hist(RitpRCw(:,2));
                title('RitpCw');

                waitforbuttonpress();
            end


        end
        end

    end
    function [Ximg]=disparity_contrast_by_CPs(cpLookup,IszRC,kernSz,k,winWHm,winPosXYZm,vDisp,varargin)
        %DB:    LExyz,RExyz,IszRC,cpLookup,
        %VDisp: IppXm,IppZm,Xpix,Yix
        %Param: PszRC,k,

        kernSz=kernSz+double(mod(kernSz,2)==0);
        S=struct(varargin{:});
        if isfield(S,'dnk') && S.dnk ~=1
            kernSzDnk=round(kernSz/S.dnk);
            kernSzDnk=kernSzDnk+double(mod(kernSzDnk,2)==0);
        else
            kernSzDnk=kernSz;
        end
        Opts=XYZ.parse_disparity_contrast__(kernSzDnk,varargin{:});
        % ENSURE ODD kernSz

        Opts.bConv=false;

        % M MESH
        W=winWHm(1)/2;
        H=winWHm(2)/2;
        [Xm,Ym]=meshgrid(linspace(-W,W,kernSzDnk(2)),linspace(H,-H,kernSzDnk(1)));
        if Opts.dnk==1
            XY=cat(3,Xm(:),Ym(:));
        else
            XY=cat(3,Xm,Ym);
        end

        % VRG VARS
        if numel(Opts.dim)==1 & Opts.dim==2
            C=0;
        else
            C=Opts.IPDm.^2;
        end

        ptchNumelDnk=prod(kernSzDnk);
        sz=[1,numel(Opts.dim)+1];
        if Opts.dim==0
            LExyz=vDisp.SubjInfo.LExyz;
            RExyz=vDisp.SubjInfo.RExyz;
        else
            LExyz=reshape(vDisp.SubjInfo.LExyz([Opts.dim 3]),sz);
            RExyz=reshape(vDisp.SubjInfo.RExyz([Opts.dim 3]),sz);
            LExyz=repmat(LExyz,ptchNumelDnk,1);
            RExyz=repmat(RExyz,ptchNumelDnk,1);
        end
        vrgCtr=2*atand(Opts.IPDm/(2*winPosXYZm(3))); % 3.7229 deg for wheatstone

        % FORWARD INTERP
        ptchNumel=prod(kernSz);
        LExyzVec=repmat(vDisp.SubjInfo.LExyz,ptchNumel,1);
        RExyzVec=repmat(vDisp.SubjInfo.RExyz,ptchNumel,1);
        Zvec=repmat(vDisp.Zm,ptchNumel,1);
        cinit=zeros(size(RExyzVec));

        % WINDOW
        WW=Opts.W(:);
        wSum=sum(WW);


        % ITER
        Ximg=nan(IszRC);
        ctrRC=floor(kernSz./2+1);
        ctr=sub2ind(kernSz,ctrRC(1),ctrRC(2));
        Rh=floor(kernSz(1)/2);
        Rw=floor(kernSz(2)/2);
        ROWS = Rh+1:IszRC(1)-Rh-1;
        COLS = Rw+1:IszRC(2)-Rw-1;

        % Prog
        N=numel(ROWS)*numel(COLS);
        modIter=floor(N/100);
        %pr=Pr(N,modIter,'disparity contrast img');
        % 5 min per patch -> 7 hours for set


        % WINDOW M TO PIX
        % XXX check XYZ back/forward project
        WszRCpix=fliplr(winWHm.*vDisp.pixPerMxy);
        [winPosRCpixL,winPosRCpixR] =vDisp.PP.back_project(winPosXYZm,true);

        AitpRC0=CPs.getAitpRC0(kernSz);


        % PATCH TO WIN SCALE
        %cpLookup=cpLookup{k};

        for i = ROWS
        for j = COLS
            ActrRC=[i,j];


            % scene to floating patch
            [LitpRC,RitpRC,LctrRC,RctrRC]=CPs.getPatchFast(ActrRC,AitpRC0,k,cpLookup,IszRC,true); % XXX Bottleneck 2 (12%)
            if isempty(LitpRC) || any(isnan(LitpRC),'all') || any(isnan(RitpRC),'all')
                continue
            end

            % float patch to centered disp.win
            [LitpRCw,RitpRCw]=CPs.patchPixToWinPix(LitpRC,RitpRC,WszRCpix,kernSz,winPosRCpixL,winPosRCpixR);

            % centered win to XYZ
            xyz=vDisp.PP.forward_project(LitpRCw,RitpRCw, LExyzVec,RExyzVec,Zvec,cinit); % XXX Bottleneck 1 (74%) : interp, intersectLinesFromPoints


            % DOWNSAMPLE HERE
            if Opts.dnk~=1
                xyz=reshape(xyz,[kernSz 3]);
                z=imresize(xyz(:,:,3),kernSzDnk,'bilinear');
                xyz=cat(3,XY,z);

                vrgImg = XYZ.to_vrg_vec_(C,LExyz,RExyz,xyz,Opts.dim);
                dsp=vrgImg(:)-vrgCtr;
            else
                vrgImg = XYZ.list_to_vrg_vec_(C,LExyz,RExyz,xyz,Opts.dim);
                dsp=vrgImg-vrgCtr;
                %xyz=cat(3,XY,xyz(:,3));
            end

            DC   = sum(dsp.*WW)./wSum;
            Ximg(i,j)  = single(real(sqrt( sum(WW.*( dsp-DC ).^2)./wSum )));


            bPlot=false; %NOTE
            if bPlot
                disp(['ActrRC ' Num.toStr(ActrRC) newline ...
                      'LctrRC ' Num.toStr(LctrRC) newline ...
                      'RctrRC ' Num.toStr(RctrRC) newline ]);

                Fig.newU(1);
                fr=4;

                subplot(fr,3,1)
                Plot.RC3(xyz,'.k');
                colormap hot;
                title('XYZm');
                ylim([.85 1.15]);
                axis square;

                subplot(fr,3,2)
                imagesc(reshape(vrgImg,Opts.PszRC)');
                set(gca,'YDir','normal');
                Fig.formatIm;
                colormap hot;
                colorbar;
                Fig.format('','','vrg deg');

                subplot(fr,3,3)
                imagesc(reshape(dsp,Opts.PszRC)');
                set(gca,'YDir','normal');
                colorbar;
                Fig.formatIm;
                colormap hot;
                Fig.format('','','dsp deg');

                % ---

                subplot(fr,3,4)
                hist(xyz(:,1));
                title('x');

                subplot(fr,3,5)
                hist(xyz(:,2));
                title('y');

                subplot(fr,3,6)
                hist(xyz(:,3));
                title('z');

                % ---

                subplot(fr,3,7)
                hist(LitpRC(:,1));
                title('LitpR');

                subplot(fr,3,8)
                hist(RitpRC(:,1));
                title('RitpR');

                subplot(fr,3,9)
                hist(RitpRC(:,2));
                title('RitpC');


                % ---

                subplot(fr,3,10)
                hist(LitpRCw(:,1));
                title('LitpRw');

                subplot(fr,3,11)
                hist(RitpRCw(:,1));
                title('RitpRw');

                subplot(fr,3,12)
                hist(RitpRCw(:,2));
                title('RitpCw');

                waitforbuttonpress();
            end


        end
        end

    end
    function [Ximg,Ximg2]=disparity_contrast_by_xyz(xyz,kernSz,winPosXYZm,winWHm,varargin)
        %method1 downsample original xyz img
        %method2 downsample xyz patch
        method=1;


        %% OPTS
        %%   dnk
        %%   dim
        %%   IPDdm
        %%   LorRorC
        %%   W
        %%   Wk
        PszRCD=size(xyz);

        % METHOD 1
        if method==1;
            S=struct(varargin{:});
            if isfield(S,'dnk') && S.dnk ~=1
                xyz=imresize(xyz,1/S.dnk,'bilinear');
                kernSz=round(kernSz/S.dnk);
            end
        end
        kernSz=kernSz+double(mod(kernSz,2)==0);

        Opts=XYZ.parse_disparity_contrast__(kernSz,varargin{:});
        Opts.bConv=false;


        ctrRC=floor(kernSz./2+1);
        ctr=sub2ind(kernSz,ctrRC(1),ctrRC(2));
        Rh=floor(kernSz(1)/2);
        Rw=floor(kernSz(2)/2);


        i=Rh+1;
        j=Rw+1;
        I=xyz(Rw+1-Rh:(Rw+1+Rh), j-Rw:j+Rw,:);

        W=winWHm(1)/2;
        H=winWHm(2)/2;
        [X,Y]=meshgrid(linspace(-W,W,kernSz(2)),linspace(H,-H,kernSz(1)));

        % METHOD 2
        if method==2 && Opts.dnk~=1
            Opts.W=imresize(Opts.W,1/Opts.dnk,'bilinear');
            I=imresize(I(:,:,1),1/Opts.dnk,'bilinear');
            I=cat(3,I,I,I);
            X=imresize(X,1/Opts.dnk,'bilinear');
            Y=imresize(Y,1/Opts.dnk,'bilinear');
        end
        WW=Opts.W(:);
        wSum=sum(WW);
        [C,LExyz,RExyz]=XYZ.parse_to_vrg__(Opts.LorRorC,Opts.IPDm,I,Opts.dim,true);
        C=C.^2;
        XY=cat(3,X,Y);

        sz=size(xyz);
        Ximg=nan(sz(1:2));

        ROWS = Rh+1:size(xyz,1)-Rh-1;
        COLS = Rw+1:size(xyz,2)-Rw-1;
        N=numel(ROWS)*numel(COLS);


        pCtr=floor(kernSz/2+1); % verified

        modIter=floor(N/100);
        %pr=Pr(N,modIter,'disparity contrast img');
        % 5 min per patch -> 7 hours for set
        vrgCtr=2*atand(Opts.IPDm/(2*winPosXYZm(3))); % 3.7229 deg for wheatstone
        for i = ROWS
            for j=COLS
                %pr.u();

                % CROP
                I=xyz((i-Rh):(i+Rh), (j-Rw):(j+Rw), :);

                % scale z by how much x or y scales
                mult=XYZ.get_scale_mult_z__(I,winWHm,kernSz,kernSz); %% BOTTLNECK 2 (getWHD) -> 24%

                %% SCALE & CENTER Z
                z0=I(pCtr(1),pCtr(2),3);
                z=((  I(:,:,3)  - z0)./mult) + winPosXYZm(3); % SLOW
                %z=I(:,:,3);

                % METHOD 2
                if method==2 && Opts.dnk ~=1
                    z=imresize(z,1/Opts.dnk,'bilinear');
                end

                % SHEER (, scale, and center) X and Y
                tI=cat(3,XY,z);

                % horizontal vrgImg
                vrgImg = XYZ.to_vrg_vec_(C,LExyz,RExyz,tI,Opts.dim); %% BOTTLENECK 1 (SSSdegSqr, A, B) -> 50%
                dsp=vrgImg(:)-vrgCtr;


                % disparityContrastImage
                DC   = sum(dsp.*WW)./wSum;
                Ximg(i,j)  = single(real(sqrt( sum(WW.*( dsp-DC ).^2)./wSum )));
            end
        end
        % METHOD 1
        if method==1 && Opts.dnk~=1
            Ximg=imresize(Ximg,PszRCD(1:2),'bilinear');
        end
        if nargout > 1
            Ximg2=rmsDeviationLoc(Ximg,WW);
        end
        %pr.c();

    end
    function Ximg=disparity_contrast(xyz,kernSz,varargin)

        if all(size(xyz(:,:,1)) > kernSz)
            xyz=Map.cropImgCtr(xyz,[],fliplr(kernSz));
        end
        Opts=XYZ.parse_disparity_contrast__(kernSz,varargin{:});
        Ximg=XYZ_disparity.disparity_contrast_(xyz,kernSz,Opts.W,Opts.dnk,Opts.IPDm,Opts.LorRorC,Opts.bConv,Opts.dim,Opts.winPosXYZm);
    end
    function Ximg=disparity_contrast_(xyz,kernSz,W,dnk,IPDm,LorRorC,bConv,dim,winPosXYZm);
        % KERNSZ DONE
        vrgCtr=2*atand(IPDm/(2*winPosXYZm(3))); % 3.7229 deg for wheatstone DONE
        %vrgImg = XYZ.toVrg(LorRorC,IPDm,xyz,dim);
        LExyz=[-IPDm/2 0 0 ]; % DONE
        RExyz=[+IPDm/2 0 0 ]; % DONE
        % dim DONE
        % C DONE
        % size(xyz) [32 32]

        vrgImg = XYZ.to_vrg_vec_(LorRorC,LExyz,RExyz,xyz,dim);
        dsp=vrgImg-vrgCtr;

        bPlot=false; %NOTE
        if bPlot
            N=5;
            subplot(1,N,1)
            imagesc(xyz(:,:,3));
            colorbar;
            title('xyz');
            axis square;
            subplot(1,N,2)
            imagesc(vrgImg);
            colorbar;
            title('vrg');
            axis square;
            subplot(1,N,3)
            imagesc(dsp);
            colorbar;
            title('dsp');
            axis square;
            subplot(1,N,4)
            imagesc(W);
            axis square;
            colorbar;

            subplot(1,N,1)
            Plot.XYZ(xyz,'.k');
            axis square;
        end
        %imagesc(W)

        if bConv==1
            Ximg   = single(real(rmsDeviationLoc(dsp,W)));
        elseif bConv==-1
            dsp=dsp(:);
            WW=W(:);

            wSum=sum(WW);
            DC   = sum(dsp.*WW)./wSum;
            Ximg = single(real(sqrt( sum(WW.*( dsp-DC ).^2)./wSum )));
            %Ximg  = single(real(sqrt( sum( W(:).*( dsp(:) ).^2 )./sum(W(:)) )));
        else
            Ximg   = single(real(rmsDeviation(dsp,W)));
        end
    end
    function Ximg=disparity_contrast_win_(xyz,kernSz,winPosXYZm,winWHm,W,dnk,IPDm,LorRorC,bConv,dim)
        xyz=XYZ.transform_to_win(xyz,winPosXYZm,winWHm,kernSz);
        Ximg=XYZ.disparity_contrast_ctr_(xyz,kernSz,W,dnk,IPDm,LorRorC,bConv,dim);
    end
    function out=disparity_contrast_ctr_(xyz,kernSz,W,dnk,IPDm,LorRorC,bConv,dim)
        vrgImg = XYZ.toVrg(LorRorC,IPDm,xyz,dim);
        out  = single(real(sqrt( sum( W(:).*( vrgImg(:) ).^2 )./sum(W(:)) )));
    end
    function Ximg=disparity_contrast_2(xyz,Wk,kernSz,IPDm,LorR,dnk)
        % XXX NEEDS WORK
        % NOTE in arcmin
        %
        W      = cosWindow(kernSz,Wk/100); %W100
        W      = W./sum(W(:));
        vrgImg = 60*vergenceFromRangeXYZVec(LorR,IPDm,xyz);
        Ximg   = single(real(rmsDeviationLoc(rmsDeviationLoc(vrgImg,W),W)));
    end
%% PARSE
    function Opts=parse_disparity_contrast__(kernSz,varargin)
        if isempty(varargin)
            Opts=Struct();
        elseif isstruct(varargin{1})
            Opts=varargin{1};
        elseif Args.isPairs(varargin{:})
            Opts=struct(varargin{:});
        end
        P={ 'PszRC', [], 'Num.is'...
           ;'Wk', 0,'Num.isInt'...
           ;'W', [],'Num.is'...
           ;'dnk',1,'Num.is'...
           ;'IPDm',0.065,'Num.is'...
           ;'LorRorC','C','Str.Alph.isLorRorC'...
           ;'bConv',[],'Num.isBinary_e'...
           ;'dim',1,'Num.is'...
           ;'winPosXYZm',[],'Num.is'
          };
        Opts=Args.parse([],P,Opts);

        if isempty(Opts.PszRC)
            Opts.PszRC=kernSz;
        end
        if isempty(Opts.bConv) || ~Opts.bConv
            Opts.bConv=isequal(kernSz,Opts.PszRC)*-1;
        end

        if isempty(Opts.W) && Opts.Wk==0
            Opts.W      = ones(kernSz);
            Opts.W      = Opts.W./sum(Opts.W(:));
        elseif isempty(Opts.W)
            Opts.W      = cosWindow(kernSz,Opts.Wk/100); %W100
            Opts.W      = Opts.W./sum(Opts.W(:));
        end

        if Opts.bConv==-1
            k=floor((Opts.PszRC-kernSz)./2);
            if ~all(k==0)
                Opts.W=padarray(Opts.W,k,0,'both');
            end
        end

    end
end
end
