classdef CP_verify < handle
methods(Static)
    function [Ximg]=dot(Lpht,Rpht,kernSz,bPixelCtrs,cpLookup,IszRC,k,vDisp)
        if k==1
            Apht=Lpht;
            Bpht=Rpht;
        elseif k==2
            Apht=Rpht;
            Bpht=Lpht;
        end

        Rh=floor(kernSz(1)/2);
        Rw=floor(kernSz(2)/2);
        ROWS = Rh+1:IszRC(1)-Rh-1;
        COLS = Rw+1:IszRC(2)-Rw-1;

        AitpRC0=CPs.getAitpRC0(kernSz);

        Ximg=nan(IszRC);
        %AVal=nan(IszRC);
        %BVal=nan(IszRC);
        Bind=nan(prod(IszRC),2);
        c=0;
        PszXY=fliplr(kernSz);

        %[Y, X] = ndgrid(1:size(Apht,1),1:size(Apht,2));
        %GA=griddedInterpolant(Y,X,Apht,'linear');
        %GB=griddedInterpolant(Y,X,Bpht,'linear');
        %
        if k==1
            nk=2;
        else
            nk=1;
        end

        if bPixelCtrs
            ActrRCAll=(1:prod(IszRC))';
        else
            ActrRCAll=cpLookup{k}{k}(1:prod(IszRC),:);
        end
        BctrRCAll=cpLookup{k}{nk}(1:prod(IszRC),:);

        for i = ROWS
        for j = COLS
            c=c+1;
            ActrRC=ActrRCAll(c,:);
            BctrRC=BctrRCAll(c,:);

            AitpRC=AitpRC0+ActrRC;
            BitpRC=AitpRC0+BctrRC;
            if any(any(isnan(AitpRC) | AitpRC < 1))
                continue
            end
            if any(any(isnan(BitpRC) | BitpRC < 1))
                continue
            end

            %AVal=GA(AitpRC(:,1),AitpRC(:,2));
            AVal=Map.cropImgCtr(Apht,[],PszXY);
            if any(isnan(AVal))
                continue
            end

            %BVal=GB(BitpRC(:,1),BitpRC(:,2));
            BVal=Map.cropImgCtr(Bpht,[],PszXY);
            if any(any(isnan(BVal)))
                continue
            end

            %Ximg(c)=(AVal'*BVal(:))/(norm(AVal)*norm(BVal(:)));
            Ximg(c)=(AVal(:)'*BVal(:))/(norm(AVal(:))*norm(BVal(:)));
        end
        end
        %sum(sum(Ximg > 1))

        %nanInd=any(~isnan(Bind),2);
        %BVal(~nanInd)=interp2(X,Y,Apht,Bind(~nanInd,2),Bind(~nanInd,1));

        %Ximg(~nanInd)=arrayfun(@(a,b) dot(a,b),AVal(~nanInd),BVal(~nanInd));

    end
    function dots=pht_dot_ptchs(P,kernSz,bSave)
        P.get_patch(1);
        PszRC=P.ptch.PszRCbuff();
        [X, Y] = meshgrid(1:PszRC(2),1:PszRC(1));
        ctr=PszRC/2;
        if mod(PszRC,2)==0
            ctr=ctr+.5;
        end
        Q=CPs.getAitpRC0(kernSz)+ctr;
        QY=Q(:,1);
        QX=Q(:,2);
        dots=zeros(length(P.fnames),1);
        %pr=Pr(length(P.fnames),1);
        kernSzXY=fliplr(kernSz);
        for i = 1:length(P.fnames)
            try
                P.get_patch(i,[],true);
                %P.get_patch(i,[]rue);
            catch
                continue
            end
            Lpht=P.ptch.mapsBuff.pht{1};
            Rpht=P.ptch.mapsBuff.pht{2};

            %RVal=interp2(X,Y,Lpht,QY,QX,'linear');
            %LVal=interp2(X,Y,Rpht,QY,QX,'linear');
            LVal=Map.cropImgCtr(Lpht,[],kernSzXY);
            RVal=Map.cropImgCtr(Rpht,[],kernSzXY);
            if isequal(kernSz,[1 1])
                dots(i)=LVal/RVal;
            else
                dots(i)=(RVal(:)'*LVal(:))/(norm(RVal(:))*norm(LVal(:)));
            end
            %pr.u();
        end
        if bSave
            fname=[P.get_dir() '_CP_verify_dot_' strrep(Num.toStr(kernSz),',','-') '_'];
            save(fname,'dots');
        end
    end
    function [gdDots,bdDots,edgesD,ctrsD]=plot_dot(P,kernSz,bGd,bBd)
        % ~ .996
        fname=[P.get_dir() '_CP_verify_dot_' strrep(Num.toStr(kernSz),',','-') '_.mat'];
        S=load(fname);
        if isequal(kernSz,[1 1])
            dots=S.dots;
            gdDots=(dots(bGd));
            bdDots=(dots(bBd));
            %binsD=linspace(0,.67,50);
            binsD=linspace(-2,2,50);
            %binsD=linspace(0,1,50);
        else
            dots=-1*S.dots+1;
            %binsD=linspace(.99,1,50);
            binsD=linspace(.9,1,50);
            %binsD=50;
            gdDots=(dots(bGd));
            bdDots=(dots(bBd));
            binsD=fliplr(-1*binsD)+1;
        end
        [~,ctrsD]=hist(gdDots,binsD);

        subPlot([1 2],1,1);
        hist(gdDots,binsD);
        if numel(binsD) > 1
            xlim([binsD(1) binsD(end)]);
        end
        title('bGd');

        subPlot([1 2],1,2);
        hist(bdDots,binsD);
        if numel(binsD) > 1
            xlim([binsD(1) binsD(end)]);
        end
        title('bBd');
        edgesD=Hist.ctrs2edges(ctrsD);
    end
    function [gdRms,bdRms,edgesR,ctrsR]=plot_rms(P,bGd,bBd)
        % ~ .06
        binsR=50;
        binsR=linspace(0.06,.2,50);
        binsR=linspace(0,.2,50);

        % RMS
        fname=[P.get_dir() '_RMSmono_.mat'];
        S=load(fname);
        rmsMonoIm=S.rmsMonoIm;

        rms=max(rmsMonoIm,[],2);

        %rms=mean(rmsMonoIm,2);

        gdRms=rms(bGd);
        bdRms=rms(bBd);


        [~,ctrsR]=hist(gdRms,binsR);
        edgesR=Hist.ctrs2edges(ctrsR);

        subPlot([1 2],1,1);
        hist(gdRms,binsR);
        title('Gd Im');
        if numel(binsR) > 1
            xlim([binsR(1) binsR(end)]);
        end

        subPlot([1 2],1,2);
        hist(bdRms,binsR);
        title('Bad Im');
        if numel(binsR) > 1
            xlim([binsR(1) binsR(end)]);
        end

    end
    function [dvDot,dvRms]=view_dot_ptchs(P,kernSz,bins)
        if nargin < 1
            P=ptchs.getRaw('all');
        end

        [bGd,bBd]=CP_verify.get_GdBd(P);
        if nargin >=3 && ~isempty(bins)
            bBin=ismember(P.idx.B,bins);
            bGd=bGd & bBin;
            bBd=bBd & bBin;
        end

        % DOT
        figure(1)
        [gd1,bd1,edges1,ctrs1]=CP_verify.plot_dot(P,kernSz,bGd,bBd);

        figure(2)
        [gd2,bd2,edges2,ctrs2]=CP_verify.plot_rms(P,bGd,bBd);

        %figure(1)
        %[gd2,bd2,edges2]=CP_verify.plot_dot(P,[3,3],bGd,bBd);
        %binsD=50;


        [N,E1,E2,loc]=histcounts2(gd1,gd2,edges1,edges2);
        figure(3)
        subPlot([1 2],1,1);
        ImapPrb.plot_fun(N,E1,E2,false);
        title('Good');

        [N,E1,E2,loc]=histcounts2(bd1,bd2,edges1,edges2);
        subPlot([1 2],1,2);
        ImapPrb.plot_fun(N,E1,E2,false);
        title('Bad');

        x=[];
        fit1='pareto';
        mDotGd=measDist('dotGd',edges1,gd1,ctrs1,fit1);
        mDotBd=measDist('dotBd',edges1,bd1,ctrs1,fit1);
        dvDot=dvDists(mDotGd,mDotBd);
        dvDot.fignum=400;
        dvDot.plot;
        dvDot.get_IN

        fit2='weib';
        mRmsGd=measDist('rmsGd',edges2,gd2,ctrs2,fit2);
        mRmsBd=measDist('rmsbD',edges2,bd2,ctrs2,fit2);
        dvRms=dvDists(mRmsGd,mRmsBd);
        dvRms.fignum=500;
        dvRms.plot;

    end
    function rms_ptchs(P,bSave)
        if nargin < 2 || isempty(bSave)
            bSave=true;
        end
        if nargin < 1
            P=ptchs.getBlk('raw');
        end
        [bGd,bBd]=CP_verify.get_GdBd(P);
        Pinds=find(bGd | bBd);

        rmsMono=nan(size(P.fnames,1),2);
        rmsMonoIm=nan(size(P.fnames,1),2);
        P.get_patch(1);
        crpSz=ceil(P.ptch.PszRC/3);
        W=ones(crpSz);
        opts=P.ptchOpts;
        opts.trgtInfo.trgtDsp=0;
        for i = 1:length(Pinds)
            try
                P.get_patch(Pinds(i),[],true);
            catch
                continue
            end
            p=P.apply_ptchOpts(P.ptch,Pinds(i),opts);
            rmsMono(Pinds(i),:)=p.im.RMSmono;
            imL=Map.cropImgCtr(p.im.img{1},[],crpSz);
            imR=Map.cropImgCtr(p.im.img{2},[],crpSz);
            rmsMonoIm(Pinds(i),1)=Map.rmsContrast(imL,W);
            rmsMonoIm(Pinds(i),2)=Map.rmsContrast(imR,W);
        end

        OUT=struct;
        OUT.rmsMono=rmsMono;
        OUT.rmsMonoIm=rmsMonoIm;
        if bSave
            fname=[P.get_dir() '_RMSmono_.mat'];
            save(fname,'rmsMono','rmsMonoIm');
        end
    end
    function [bGd,bBd]=get_GdBd(P,bin)
        if nargin < 2 || isempty(bin)
            bInd=true(size(P.fnames));
        else
            bInd=P.idx.B==bin;
        end
        nInd=main.newInd(P);

        bBd=~nInd & bInd &(P.Flags.seen | P.Flags.other) &  P.Flags.bad;
        bGd=~nInd & bInd &(P.Flags.seen | P.Flags.other) & ~P.Flags.bad;
    end
    function view_rms_ptchs(P)

    end
end
end
