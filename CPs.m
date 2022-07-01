classdef CPs < handle & CPs_exp
properties(Constant)
    LANDR='LRC'
    NOTLANDR='RL'
    NOTK=[2 1]
end
methods(Static)
%%- UTIL
    function [k,nk]=getK(LorRorK)
        if isnumeric(LorRorK)
            k=LorRorK;
        elseif LorRorK=='L'
            k=1;
        elseif LorRorK=='R'
            k=2;
        end

        if k==1
            nk=2;
        elseif k==2
            nk=1;
        end
    end
 %% THIS
    function [LorR,notLorR]=getLorR(k)
        if ischar(k)
            LorR=k;
        elseif k==3
            LorR='C';
        else
            LorR=char(transpose(CPs.LANDR(k)));
        end
        notLorR=char(transpose(CPs.NOTLANDR(k)));
    end
    function ind=getOutOfRangeInd(ActrRC,IszRC)
        ind=any(isnan(ActrRC) | ActrRC<=0,2) | ActrRC(:,1) > IszRC(1) | ActrRC(:,2) > IszRC(2);
    end
    function ind=getIndOutOfRangeInd(ActrInd,IszRC)
        ind=isnan(ActrInd) | ActrInd > prod(IszRC) | ActrInd <= 0;
    end
    function ind=nearest(ptRC,IszRC)
        r=round(ptRC);
        r=arrayfun(@(x) max(x,1),r);
        r(:,1)=arrayfun(@(x) min(x,IszRC(1)),r(:,1));
        r(:,2)=arrayfun(@(x) min(x,IszRC(2)),r(:,2));
    end
    function idealRC=getAllRC(IszRC)
        idealRC=fliplr(Set.distribute(1:IszRC(2),1:IszRC(1)));
    end
    function ind=getAllInd(IszRC)
        ind=transpose(1:prod(IszRC));
    end %THIS IS
    % this
    %
    function [xOffset,d]=getXOffset(AitpRC,IszRC)
        IitpRC=CPs.getAllRC(IszRC);
        d=AitpRC(:,2)-IitpRC(:,2);
        xOffset=nanmedian(d);
    end
    function AitpRC0=getAitpRC0(PszRC);
        [x0,y0]=Rec.rectPix([0,0],PszRC(1),PszRC(2));
        AitpRC0=Set.distribute(y0(3):y0(1), x0(2):x0(1));
    end
%%- FNAME
    function name=getName(I)
        name=[num2str(I,'%03i')]; %
        name=[num2str(I,'%03i')];
    end
    function fname=getFname(database,I)
        dire=CPs.getDire(database);
        name=CPs.getName(I);
        fname=[dire name];
    end
    function dire=getDire(database)
        dire=Env.var('ImgDb.CPs',database);
    end
    function cps=load(database,I,k)
        fname=CPs.getFname(database,I);
        S=load(fname);
        cps=S.CPs;

    end
%% FULL
    function name=getFullName(I,k)
        if ischar(i)
            LorR=k;
        elseif Num.is(k)
            LorR=CPs.LANDR(k);;
        end
        name=[LorR num2str(I,'%03i')];
    end
    function dire=getFullDire(database)
        dire=Env.var('ImgDb.CPsFull',database);
    end
    function fname=getFullFname(database,I,K)
        dire=CPs.getFullDire(database);
        name=CPs.getFullName(I,K);
        fname=[dire name];
    end
    function cps=loadFull(database,I,k, bBaseLoad)
        if nargin < 4
            bBaseLoad=false;
        end
        fname=CPs.getFullFname(database,I,k);
        if bBaseLoad
            k=CPs.getK(k);
            varName=['cps_' num2str(I) '_' num2str(k)];
            Data.baseLoad(fname,varName);
        else
            load(fname);
        end
    end
%% LOOKUP
    function name=getLookupName(I,bChk)
        name=[num2str(I,'%03i')];
        if ~exist('bChk','var') || isempty(bChk) || ~bChk
            return
        end
        name=['chk' name];
    end
    function dire=getLookupDire(database)
        dire=Env.var('ImgDb.CPsLookup',database);
    end
    function fname=getLookupFname(database,I,bChk)
        if ~exist('bChk','var') || isempty(bChk)
            bChk=false;
        end
        dire=CPs.getLookupDire(database);
        name=CPs.getLookupName(I,bChk);
        fname=[dire name '.mat'];
    end
    function cpLookup=loadLookup(database,I,bChk,bBaseLoad)
        if nargin < 3 || isempty(bChk)
            bChk=false;
        elseif ischar(bChk)
            error('CPs.oadLookup does not take LorR param');
        end
        if nargin < 4
            bBaseLoad=false;
        end
        fname=CPs.getLookupFname(database,I,bChk);
        if bBaseLoad
            varName=['cps_' num2str(I)];
            Data.baseLoad(fname,varName);
        else
            load(fname);
        end
        if bChk
            cpLookup=cpLookupChk;
        end
    end

%% XYZ
    function name=getXYZName(I,bChk)
        name=[num2str(I,'%03i')];
        if ~exist('bChk','var') || isempty(bChk) || ~bChk
            return
        end
        name=['chk' name];
    end
    function dire=getXYZDire(database)
        dire=Env.var('ImgDb.CPsXYZ',database);
    end
    function fname=getXYZFname(database,I,bChk)
        if ~exist('bChk','var') || isempty(bChk)
            bChk=false;
        end
        dire=CPs.getXYZDire(database);
        name=CPs.getXYZName(I,bChk);
        fname=[dire name '.mat'];
    end
    function cpLookup=loadXYZ(database,I,bChk)
        if nargin < 3 || isempty(bChk)
            bChk=false;
        end

        fname=CPs.getXYZFname(database,I,bChk);
        load(fname);
    end
%%- LOOKUP
    function [LitpRC,RitpRC]=lookupLR(AitpRC,k,cpLookup,IszRC, bFlip)
        bFlip=0;
        if nargin < 5
            bFlip=false;
        end
        [k,nk]=CPs.getK(k);

        AitpRC=CPs.roundcheck(AitpRC,IszRC);
        if size(AitpRC,2) == 2
            inds=Sub.toInd(AitpRC,IszRC);
            %inds=sub2ind(IszRC,AitpRC(:,1),AitpRC(:,2));
        elseif size(AitpRC,2) == 1
            inds=AitpRC;
        end
        %size(cpLookup{k}{1}(inds(1),:))
        %size(cpLookup{k}{2})

        %k

        LitpRC=cpLookup{k}{1}(inds,:);
        RitpRC=cpLookup{k}{2}(inds,:);

        % XXX
        %if bFlip
        %    LitpRC=fliplr(LitpRC);
        %    RitpRC=fliplr(RitpRC);
        %end
        %
        if bFlip
            if k==1
                [LitpRC,RitpRC]=CPs.flipcheck(LitpRC,RitpRC);
            else
                [RitpRC,LitpRC]=CPs.flipcheck(RitpRC,LitpRC);
            end
        end

        %[max(AitpRC(:,1)) max(AitpRC(:,2))]
        %if k == 1
        %    [max(RitpRC(:,1)) max(RitpRC(:,2))]
        %else
        %    [max(LitpRC(:,1)) max(LitpRC(:,2))]
        %end

    end
    function [AitpRC,BitpRC]=lookupAB(AitpRC,k,cpLookup,IszRC,bFlip)
        bFlip=0;
        if nargin < 5
            bFlip=false;
        end
        [k,nk]=CPs.getK(k);

        AitpRC=CPs.roundcheck(AitpRC,IszRC);

        inds=Sub.toInd(AitpRC,IszRC);
        AitpRC=cpLookup{k}{k}(inds,:);
        BitpRC=cpLookup{k}{nk}(inds,:);

        % % XXX
        % if abs(mod(AitpRC(1,1),1))>0
        %     if size(AitpRC,1) > 1
        %         AitpRC(5,:)
        %     end
        %     AitpRC=fliplr(AitpRC);
        %     if size(AitpRC,1) > 1
        %         AitpRC(5,:)
        %     end
        % end

        if bFlip
            [AitpRC,BitpRC]=CPs.flipcheck(AitpRC,BitpRC);
        end
    end
    function AitpRC=roundcheck(AitpRC,IszRC)
        AitpRC=round(AitpRC); %% NOTE
        lrg=AitpRC == IszRC;
        sml=AitpRC == 0;
        if any(lrg(:,1))
            AitpRC(lrg(:,1))=IszRC(1);
        end
        if any(lrg(:,2))
            AitpRC(lrg(:,2))=IszRC(2);
        end
        if any(sml)
            AitpRC(sml)=1;
        end

    end
    function [AitpRC,BitpRC]=flipcheck(AitpRC,BitpRC)
        if abs(mod(BitpRC(1,1),1))>0
             BitpRC=fliplr(BitpRC);
        end
        if abs(mod(AitpRC(1,1),1))>0 && (abs(BitpRC(1,1) - AitpRC(1,1)) > 2) && (abs(BitpRC(1,1) - AitpRC(1,2)) < 2)
            AitpRC=fliplr(AitpRC);
        end
    end
    function [AitpRC,BitpRC]=lookupQ(AitpRC,k,cpLookup,IszRC)
        [k,nk]=CPs.getK(k);
        inds=Sub.toInd(IszRC,AitpRC);
        AitpRC=cpLookup{k}{1}(inds,:);
        BitpRC=cpLookup{nk}{2}(inds,:);
    end
    function [AitpInd,BitpInd,LctrRC,RctrRC]=ActrRCToABRecInd(ActrRC,AitpRC0,k,cpLookup,IszRC);
        [LctrRC,RctrRC]=CPs.lookupLR(ActrRC,k,cpLookup,IszRC);

        if k==1
            AitpRC=AitpRC0 + LctrRC;
            BitpRC=AitpRC0 + RctrRC;
        else
            AitpRC=AitpRC0 + RctrRC;
            BitpRC=AitpRC0 + LctrRC;
        end
        AitpRC=round(AitpRC);
        BitpRC=round(BitpRC);

        AitpInd=Sub.toInd(AitpRC,IszRC);
        BitpInd=Sub.toInd(BitpRC,IszRC);
    end
%%-
    function pointsXYZ=fullToXYZ(database,I,LorR,bPlot,bDebug)
        if ~exist('bPlot','var') || isempty(bPlot)
            bPlot=false;
        end
        db=dbInfo(database);
        k=CPs.getK(LorR);
        K=k;
        db=dbInfo(database);
        full=CPs.loadFull(database,I,k,bDebug);
        [X,Y]=meshgrid(1:db.IszRC(2),1:db.IszRC(1));
        pointsXYZ=XYZ.forward_project(db.LExyz,db.RExyz,fliplr(full.LitpRC),fliplr(full.RitpRC),db.IppXm{K},db.IppYm{K},db.IppZm,X,Y);
        %pointsXYZ=XYZ.forward_project(db.LExyz,db.RExyz,fliplr(LitpRCchk),fliplr(RitpRCchk),db.IppXm{K},db.IppYm{K},db.IppZm,X,Y);
        %
        if bPlot
            edges1= 1.0e+03 .* [0.0977 ,0.2910 ,0.4843 ,0.6776 ,0.8709 ,1.0642 ,1.2575 ,1.4508 ,1.6442 ,1.8375];
            edges2 = 1.0e+03 .* [0.0969 ,0.2886 ,0.4803 ,0.6721 ,0.8638 ,1.0556 ,1.2473 ,1.4390 ,1.6308 ,1.8225];

            %Fig.newU(1);
            %hist(full.LitpRC,edges1);

            %Fig.newU(2);
            %hist(full.RitpRC,edges2);

            Fig.newU(3);
            Plot.RC3(pointsXYZ,'k.');


            %Fig.newU(4);
            %Plot.imagescRC3(pointsXYZ);

            %xlim([-1 1]);
            %zlim([-0.5 0.5]);
            %ylim([1 3]);
            %imagesc(img);
            %Fig.formatIm();
            %Plot.RC3(pointsXYZ,'k.');
        end
    end

    function pointsXYZ=toXYZ(database,I,LorR,bPlot,bDebug)
        if ~exist('bPlot','var') || isempty(bPlot)
            bPlot=false;
        end
        db=dbInfo(database);
        [k,nk]=CPs.getK(LorR);
        K=3;
        db=dbInfo(database);
        lookup=CPs.loadXYZ(database,I,0,1); % XXX
        [X,Y]=meshgrid(1:db.IszRC(2),1:db.IszRC(1));

        AitpRC=CPs.getAllRC(db.IszRC);
        [LitpRC,RitpRC]=CPs.lookupAB(AitpRC,k,lookup,db.IszRC);
        %if k==1
        %    LitpRC=AitpRC;
        %elseif k==2
        %    RitpRC=AipRC;
        %end

        pointsXYZ=XYZ.forward_project(db.LExyz,db.RExyz,fliplr(LitpRC),fliplr(RitpRC),db.IppXm{K},db.IppYm{K},db.IppZm,X,Y);
        %if k ==1
        %img=F(db.IppXm{K},db.IppYm{K});
        if bPlot
            edges1= 1.0e+03 .* [0.0977 ,0.2910 ,0.4843 ,0.6776 ,0.8709 ,1.0642 ,1.2575 ,1.4508 ,1.6442 ,1.8375];
            edges2 = 1.0e+03 .* [0.0969 ,0.2886 ,0.4803 ,0.6721 ,0.8638 ,1.0556 ,1.2473 ,1.4390 ,1.6308 ,1.8225];

            %Fig.newU(1);
            %hist(lookup{2},edges1);

            %Fig.newU(2);
            %hist(lookup{1},edges2);

            Fig.newU(3);
            Plot.RC3(pointsXYZ,'k.');
            %xlim([-2 2]);

            %Fig.newU(4);
            %Plot.imagescRC3(pointsXYZ);

        end
    end
    function pointsXYZ=lookupToXYZ(database,I,LorR,bPlot,bDebug)
        if ~exist('bPlot','var') || isempty(bPlot)
            bPlot=false;
        end
        db=dbInfo(database);
        [k,nk]=CPs.getK(LorR);
        K=k;
        db=dbInfo(database);
        lookup=CPs.loadXYZ(database,I,false,bDebug); % XXX
        [X,Y]=meshgrid(1:db.IszRC(2),1:db.IszRC(1));

        if k==1
            LitpRC=CPs.getAllRC(db.IszRC);
            [LitpRC,RitpRC]=CPs.lookupLR(LitpRC,k,lookup,db.IszRC);
        elseif k==2
            RitpRC=CPs.getAllRC(db.IszRC);
            [LitpRC,RitpRC]=CPs.lookupLR(RitpRC,k,lookup,db.IszRC);
        end

        pointsXYZ=XYZ.forward_project(db.LExyz,db.RExyz,fliplr(LitpRC),fliplr(RitpRC),db.IppXm{K},db.IppYm{K},db.IppZm,X,Y);
        %if k ==1
        %img=F(db.IppXm{K},db.IppYm{K});
        if bPlot
            edges1= 1.0e+03 .* [0.0977 ,0.2910 ,0.4843 ,0.6776 ,0.8709 ,1.0642 ,1.2575 ,1.4508 ,1.6442 ,1.8375];
            edges2 = 1.0e+03 .* [0.0969 ,0.2886 ,0.4803 ,0.6721 ,0.8638 ,1.0556 ,1.2473 ,1.4390 ,1.6308 ,1.8225];

            %Fig.newU(1);
            %hist(lookup{2},edges1);

            %Fig.newU(2);
            %hist(lookup{1},edges2);

            Fig.newU(3);
            Plot.RC3(pointsXYZ,'k.');

            %Fig.newU(4);
            %Plot.imagescRC3(pointsXYZ);

        end
    end
%%- PATCH
    function [LitpRC,RitpRC]=getPatch(database,I,LorR,ActrRC,PszRC,bCtr,bChk)
        %[ LitpRC,RitpRC,AitpRC,BitpRC,BctrRC]=getPatch(database,I,LorR,ActrRC,PszRC,bCtr,bChk)

        db=dbInfo(database);
        k=CPs.getK(LorR);
        AitpRC0=CPs.getAitpRC0(PszRC);
        cpLookup=CPs.loadXYZ(database,I,bChk); % XXX

        [LitpRC,RitpRC]=CPs.getPatchFast(ActrRC,AitpRC0,k,cpLookup,db.IszRC,bCtr);
    end
    function [LitpRC,RitpRC,LctrRC,RctrRC]=getPatchFast(ActrRC,AitpRC0,k,cpLookup,IszRC,bCtr)
        AitpRC=AitpRC0+ActrRC;
        %if any(AitpRC < 0,'all')
        if any(any(AitpRC < 0))
            LitpRC=[]; RitpRC=[]; LctrRC=[]; RctrRC=[];
            return
        end

        [LctrRC,RctrRC]=CPs.lookupLR(ActrRC,k,cpLookup,IszRC,true);

        % LctrRC
        % RctrRC
        % ActrRC

        % inds=Sub.toInd(ActrRC,IszRC);
        % ActrRC
        % inds2=sub2ind(IszRC,ActrRC(1),ActrRC(2));
        % [a,b]=ind2sub(IszRC,inds);
        % [c,d]=ind2sub(IszRC,inds2);
        % subs(:,1)=a;
        % subs(:,2)=b;
        % subs2(:,1)=c;
        % subs2(:,2)=d;
        % subs
        % subs2



        % if abs(ActrRC(1)-LctrRC(2)) < 5
        %     bFlipL=true;
        %     LctrRC=fliplr(LctrRC);
        % else
        %     bFlipL=false;
        % end
        % if abs(ActrRC(1)-RctrRC(2)) < 5
        %     bFlipR=true;
        %     RctrRC=fliplr(RctrRC);
        % else
        %     bFlipR=false;
        % end


        if any(isnan(LctrRC)) || any(isnan(LctrRC))
            LitpRC=[]; RitpRC=[];
            return
        end
        [LitpRC,RitpRC]=CPs.lookupLR(AitpRC,k,cpLookup,IszRC,true);

        % if bFlipL
        %     LitpRC=fliplr(LitpRC);
        % end
        % if bFlipR
        %     RitpRC=fliplr(RitpRC);
        % end

        %33
        %[LctrRC; RctrRC]

        % subplot(1,5,1)
        % Plot.RC(AitpRC,'.');
        % subplot(1,5,2)
        % Plot.RC(LitpRC,'.');
        % subplot(1,5,3)
        % Plot.RC(RitpRC,'.');
        % size(AitpRC)
        % IszRC
        % ActrRC-LctrRC
        % ActrRC-RctrRC
        % d1=AitpRC-LitpRC;
        % d2=AitpRC-RitpRC;
        % max(d1(:,1))
        % max(d1(:,2))
        % max(d2(:,1))
        % max(d2(:,2))
        % LctrRC
        % RctrRC
        % k


        % NOTE
        %if k == 1
        %    LitpRC=AitpRC;
        %    LctrRC=ActrRC;
        %else
        %    RitpRC=AitpRC;
        %    RctrRC=ActrRC;
        %end

        if bCtr
            LitpRC=LitpRC-LctrRC;
            RitpRC=RitpRC-RctrRC;
        end

        % subplot(1,5,4)
        % Plot.RC(LitpRC,'.')
        % subplot(1,5,5)
        % Plot.RC(RitpRC,'.')
        % waitforbuttonpress

    end
    function [LitpRC,RitpRC]=patchPixToWinPix(LitpRC,RitpRC,WszRCpix,PszRCpix,WinPosRCpixL,WinPosRCpixR)
    %function [LitpRC,RitpRC]=CPs.patchPixToWinPix(LitpRC,RitpRC,WszRC,PszRC,WinPosRCpixL,WinPosRCpixR)
    % NEEDS TO BE CENTERED
        mult=WszRCpix./PszRCpix;
        % floating patch to floating disp.win to centered disp.win
        LitpRC=(LitpRC.*mult)+WinPosRCpixL;
        RitpRC=(RitpRC.*mult)+WinPosRCpixR;
    end
    function [LitpRC,RitpRC]=patchPixToWinPixBuff(LitpRC,RitpRC,WszRCpix,PszRCpix,PszRCpixBuff,WinPosRCpixL,WinPosRCpixR)

        WszRCpixBuff=WszRCpix.*PszRCpixBuff./PszRCpix;


        [LitpRC,RitpRC]=CPs.patchPixToWinPix(LitpRC,RitpRC,WszRCpixBuff,PszRCpixBuff,WinPosRCpixL,WinPosRCpixR);
    end
%%- GEN
end
methods(Static, Hidden)
%% PARAMS
%% GET
%%- UTIL
    function [LitpRCdsp,RitpRCdsp,bIndGdFXN]=addDsp(LitpRC,RitpRC,dspDeg,vDisp)
        % IppXm, center is zero
        %if ILorR=='L'
        %    LppXm=IppXm;

        %    RppXm=LppXm-IPDm;
        %elseif ILorR=='R'
        %    RppXm=IppXm;
        %    LppXm=RppXm+IPDm;
        %elseif ILorR=='C'
        %    %LppXm=IppXm+IPDm/2; % SLOW XXX
        %    %RppXm=IppXm-IPDm/2; % SLOW XXX
        %end
        %LppYm=IppXm.LppYm;
        %RppYm=IppXm.RppYm;
        [LitpRCdsp,RitpRCdsp,bIndGdFXN]=addDsp(LitpRC,RitpRC,...
                                               dspDeg,...
                                               vDisp.PP.L.Xm,...
                                               vDisp.PP.L.Ym,...
                                               vDisp.PP.R.Xm,...
                                               vDisp.PP.R.Ym,...
                                               vDisp.Zm,...
                                               vDisp.SubjInfo.IPDm);
    end
    function ppTest_()
        V=VDisp;
        IPDm=V.SubjInfo.IPDm;
        CppXm=V.PP.C.Xm;
        RppXm=V.PP.R.Xm;
        LppXm=V.PP.L.Xm;

        %imagesc((CppXm+IPDm/2)-LppXm);
        %imagesc((CppXm-IPDm/2)-RppXm);
    end
end
end
