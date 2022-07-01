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
            dots(i)=(RVal(:)'*LVal(:))/(norm(RVal(:))*norm(LVal(:)));
            %pr.u();
        end
        if bSave
            fname=[P.get_dir() '_CP_verify_dot_' strrep(Num.toStr(kernSz),',','-') '_'];
            save(fname,'dots');
        end
    end
    function dots=view_dot_ptchs(P,kernSz)
        if nargin < 1
            P=ptchs.getRaw('all');
        end
        fname=[P.get_dir() '_CP_verify_dot_' strrep(Num.toStr(kernSz),',','-') '_.mat'];
        S=load(fname);
        dots=S.dots;

        [bGd,bBd]=CP_verify.get_GdBd(P);

        gdDots=dots(bGd);
        bdDots=dots(bBd);
        bins=linspace(.99,1,50);
        %bins=linspace(0,1,50);

        figure(1)

        subPlot([1 2],1,1);
        hist(gdDots,bins);
        xlim([bins(1) bins(end)]);
        title('bGd');

        subPlot([1 2],1,2);
        hist(bdDots,bins);
        xlim([bins(1) bins(end)]);
        title('bBd');


        fname=[P.get_dir() '_RMSmono_.mat'];
        S=load(fname);
        rmsMono=S.rmsMono;
        rmsMonoIm=S.rmsMonoIm;
        mx=max(rmsMono,[],2);
        mxIm=max(rmsMonoIm,[],2);

        figure(2)
        subPlot([2 2],1,1);
        hist(mxIm(bGd));
        title('Gd Im');

        %sum(~isnan(mxIm(bBd)))
        subPlot([2 2],1,2);
        hist(mxIm(bBd));
        title('Bad Im');

        subPlot([2 2],2,1);
        hist(mx(bGd));

        subPlot([2 2],2,2);
        hist(mx(bBd));

    end
    function rms_ptchs(P,bSave)
        if nargin < 2 || isempty(bSave)
            bSave=true;
        end
        if nargin < 1
            P=ptchs.getBlk('all');
        end
        [bGd,bBd]=CP_verify.get_GdBd(P);
        Pinds=find(bGd | bBd);

        rmsMono=nan(size(P.fnames,1),2);
        rmsMonoIm=nan(size(P.fnames,1),2);
        P.get_patch(1);
        W=ones(P.ptch.PszRC);
        for i = 1:length(Pinds)
            try
                P.get_patch(i);
            catch
                continue
            end
            rmsMono(Pinds(i),:)=P.ptch.im.RMSmono;
            rmsMonoIm(Pinds(i),1)=Map.rmsContrast(P.ptch.im.img{1},W);
            rmsMonoIm(Pinds(i),2)=Map.rmsContrast(P.ptch.im.img{2},W);
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
