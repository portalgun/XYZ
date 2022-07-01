classdef CPs_db < handle
methods(Static)
    function out=xyz(database,startAt_OR_inds,bSave,bDebug)
        bAndChk=true;
        db=dbInfo(database);
        xdiff=diff(db.IppXm{1}(1,1:2));
        ydiff=abs(diff(db.IppYm{1}(1:2,1)));

        %if nargin < 5
            %method=1; %% NOTE
            if nargin < 4
                bDebug=false;
            end
            if nargin < 3
                bSave=true;
            end
        %end

        if ~exist('startAt_OR_inds','var') || isempty(startAt_OR_inds)
            INDS=1:length(db.allImages);
        elseif numel(startAt_OR_inds) == 1
            INDS=startAt_OR_inds:numel(db.allImages);

        elseif numel(startAt_OR_inds) > 1
            INDS=Vec.row(startAt_OR_inds);
            if numel(INDS)==2 && numel(unique(INDS))==1
                INDS=INDS(1);
            end
        end

        dire=CPs.getXYZDire(database);
        if ~Dir.exist(dire)
            mkdir(dire);
        end

        inds=CPs.getAllInd(db.IszRC);
        n=numel(inds);
        ActrRC=CPs.getAllRC(db.IszRC);

        Anch=3;
        [X,Y]=meshgrid(1:db.IszRC(2),1:db.IszRC(1));
        img=cell(2,1);

        N=length(INDS);
        p=Pr(N,1,'Generating CPsXYZ');
        for i = INDS
            p.u();
            I=db.allImages(i);

            cps=cell(2,1);
            cpsChk=cell(1,2);
            for k = 1:2
                [k,nk]=CPs.getK(k);
                LorR=CPs.getLorR(k);
                ii=dbImg.getImg('LRSI','img','xyz',I,k);
                %ii=permute(ii,[2 1 3]);
                img{k}=reshape(ii,[prod(db.IszRC),3]);
                full=CPs.loadFull(database,I,k,bDebug);

                if Anch==3
                    K=3;
                    LExyz=db.LExyz;
                    RExyz=db.RExyz;
                else

                    K=k;
                    if k==1
                        LExyz=db.L.AExyz;
                        RExyz=db.L.BExyz;
                    elseif k==2
                        RExyz=db.R.AExyz;
                        LExyz=db.R.BExyz;
                    end
                end

                [cps{k}{1}, cps{k}{2}]=XYZ.back_project(LExyz,RExyz,img{k},db.IppXm{K},db.IppYm{K},db.IppZm,X,Y);
                %cps{k}{1}
                %waitforbuttonpress
                %cps{k}{2}
                %waitforbuttonpress
                %dk
                %sum(mod(cps{k}{1},1)==0)
                %sum(mod(cps{k}{2},1)==0)
                %cps{k}{nk}(:,1)=cps{k}{k}(:,1);
            end

            fname=CPs.getXYZFname(database,I);
            if bDebug
                Data.baseSave(cps,['cps_' num2str(I)]);
            end
            if bSave
                cpLookup=cps;
                save(fname,'cpLookup');
            end
            if bDebug
                [k,nk]=CPs.getK(k);

                figure(1)
                subPlot([2,2],1,1);
                Plot.RC(cps{k}{1});

                subPlot([2,2],1,2);
                Plot.RC(cps{k}{2});

                subPlot([2,2],2,1);
                Plot.RC(cps{nk}{1});

                subPlot([2,2],2,2);
                Plot.RC(cps{nk}{2});

                figure(2)
                subPlot([1 2],1,1);

                if Anch==3
                    K=3;
                else
                    K=1;
                end
                pointsXYZ=XYZ.forward_project(LExyz,RExyz,cps{1}{1},cps{1}{2},db.IppXm{K},db.IppYm{K},db.IppZm,X,Y);
                Plot.RC3(pointsXYZ,'.k');

                subPlot([1 2],1,2);
                Plot.RC3(img{1},'.k');

                figure(3)
                Xerr=abs(pointsXYZ(:,2)-img{1}(:,2));
                Yerr=abs(pointsXYZ(:,1)-img{1}(:,1));
                subplot(2,1,1)
                hist(Xerr);
                subplot(2,1,2)
                hist(Yerr);

                %disp(Xerr)
                %disp(Yerr)

                drawnow
            end
        end
    end
   function out=lookup(database,startAt_OR_inds,bSave,bCheck,bDebug,method)
        bAndChk=true;
        db=dbInfo(database);
        xdiff=diff(db.IppXm{1}(1,1:2));
        ydiff=abs(diff(db.IppYm{1}(1:2,1)));

        if nargin < 6
            method=1; %% NOTE
        end
        if nargin < 5
            bDebug=false;
        end
        if nargin < 4
            bSave=true;
        end

        if ~exist('startAt_OR_inds','var') || isempty(startAt_OR_inds)
            INDS=1:length(db.allImages);
        elseif numel(startAt_OR_inds) == 1
            INDS=startAt_OR_inds:numel(db.allImages);

        elseif numel(startAt_OR_inds) > 1
            INDS=Vec.row(startAt_OR_inds);
            if numel(INDS)==2 && numel(unique(INDS))==1
                INDS=INDS(1);
            end
        end

        dire=CPs.getLookupDire(database);
        if ~Dir.exist(dire)
            mkdir(dire);
        end

        inds=CPs.getAllInd(db.IszRC);
        n=numel(inds);
        ActrRC=CPs.getAllRC(db.IszRC);

        N=length(INDS);
        p=Pr(N,1,'Generating CPsLookup');
        for i = INDS
            p.u();
            I=db.allImages(i);

            cps=cell(2,1);
            cpsChk=cell(1,2);
            for k = 1:2
                [k,nk]=CPs.getK(k);
                LorR=CPs.getLorR(k);

                full=CPs.loadFull(database,I,k,bDebug);
                LitpRC=full.LitpRC;
                RitpRC=full.RitpRC;

                cps{k}={LitpRC RitpRC};
                cpsChk{k}={full.LitpRCchk full.RitpRCchk};

            end

            if method == 2
                [cps,exitflag1,exitflag2]   =CPs_exp.toLookup(cps{1},      cps{2},db,ActrRC);
                if exitflag1 || exitflag2

                    p.append_msg(sprintf('Problem with image %d',I));
                    continue
                end
                %cpsChk=CPs_exp.toLookup(cpsChk{1},cpsChk{2},db.IszRC);
            elseif method == 3
                cps   =CPs_exp.toLookupS(cps{1},      cps{2},db.IszRC,xdiff,ydiff);
                %cpsChk=CPs_exp.toLookupS(cpsChk{1},cpsChk{2},db.IszRC);
            end
            assert(all(size(cps)==[2,1]));
            assert(all(size(cps{1})==[1,2]));
            assert(all(size(cps{2})==[1,2]));
            assert(all(size(cps{1}{1})==[n,2]));
            assert(all(size(cps{1}{2})==[n,2]));
            assert(all(size(cps{2}{1})==[n,2]));
            assert(all(size(cps{2}{2})==[n,2]));
            % NOTE
            cps{3}=full.ActrRC;


            if bDebug
                Data.baseSave(cps,['cps_' num2str(I)]);
            end

            fname=CPs.getLookupFname(database,I);
            fnameChk=CPs.getLookupFname(database,I,true);
            if bSave
                cpLookup=cps;

                save(fname,'cpLookup');

            end
            if bCheck
                AitpRC=randi(prod(db.IszRC),300,1);

                [LitpRC,RitpRC]=CPs.lookupLR(AitpRC,1,cps,db.IszRC);
                subplot(2,1,1)
                CPs_db.im_check(I,db,LitpRC,RitpRC,'all');

                [LitpRC,RitpRC]=CPs.lookupLR(AitpRC,2,cps,db.IszRC);
                subplot(2,1,2)
                CPs_db.im_check(I,db,LitpRC,RitpRC,'all');
                drawnow
            end
            if bDebug && bSave
                pointsXYZL=CPs.lookupToXYZ(database,I,1,false,bDebug);
                pointsXYZR=CPs.lookupToXYZ(database,I,2,false,bDebug);
                sum(pointsXYZL(:,3) < 0,'all')
                sum(pointsXYZR(:,3) < 0,'all')
                %find(pointsXYZL(:,3) < 0)
                %find(pointsXYZR(:,3) < 0)

                Fig.newU(1);

                subplot(2,1,1)
                Plot.RC3(pointsXYZL,'k.');
                %Plot.imagescRC3(pointsXYZL);
                axis square;

                subplot(2,1,2)
                Plot.RC3(pointsXYZR,'k.');
                %Plot.imagesc(pointsXYZL,'k.');
                axis square;

                drawnow
                waitforbuttonpress


            end
        end
        p.c();
    end
    function out=full(database,startAt_OR_inds,bSave,bCheck)
        db=dbInfo(database);
        if nargin <4
            bCheck=false;
        end
        if ~exist('startAt_OR_inds','var') || isempty(startAt_OR_inds)
            INDS=1:length(db.allImages);
        elseif numel(startAt_OR_inds) == 1
            INDS=startAt_OR_inds:numel(db.allImages);

        elseif numel(startAt_OR_inds) > 1
            INDS=Vec.row(startAt_OR_inds);
            if numel(INDS)==2 && numel(unique(INDS))==1
                INDS=INDS(1);
            end
        end

        bAndChk=true;

        ActrRC=CPs.getAllRC(db.IszRC);

        dire=CPs.getFullDire(database);
        if ~Dir.exist(dire)
            mkdir(dire);
        end

        N=length(INDS);
        p=Pr(N,1,'Generating CPsFull');
        for I =INDS
            p.u();
            [CPsL,CPsR]=CPs_gen.genImgBoth(database,I,ActrRC,bAndChk,db);
            fnameL=CPs.getFullFname(database,I,1);
            fnameR=CPs.getFullFname(database,I,2);
            if bSave
                cps=CPsL;
                save(fnameL,'cps');

                cps=CPsR;
                save(fnameR,'cps');

            end
            if bCheck
                subplot(2,1,1)
                CPs_db.im_check(I,db,CPsL.LitpRC,CPsL.RitpRC);
                subplot(2,1,2)
                CPs_db.im_check(I,db,CPsR.LitpRC,CPsR.RitpRC);
                drawnow
            end
        end
        p.c();
    end
    function im_check(I,db,LitpRC,RitpRC, inds)

        im=dbImg.getImg(db.database,'img','pht',I);
        if nargin < 5
            inds=randi(prod(db.IszRC),300,1);
        elseif ischar(inds) && strcmp(inds,'all')
            inds=1:size(LitpRC,1);
        end

        hold off
        imagesc([im{1} im{2}]);
        hold on
        Plot.RC(LitpRC(inds,:),'.r');
        Plot.RC(RitpRC(inds,:)+[0 db.IszRC(2)],'.r');
        Fig.formatIm;
    end
    function out= check_full(database,startAt)
        if ~exist('startAt','var') || isempty(startAt)
            startAt=1;
        end
        db=dbInfo(database);
        ActrRC=CPs.getAllRC(db.IszRC);
        dire=CPs.getFullDire(database);
        if ~Dir.exist(dire)
            mkdir(dire);
        end
        N=length(db.allImages);
        p=Pr(N,1,'Generating CPsFull');
        IszRC=db.IszRC;
        for i = 1:N
            p.u();
            if i < startAt
                continue
            end
            I=db.allImages(i);
            im=getImg(database,'img','pht',I);
            fnameL=CPs.getFullFname(database,I,1);
            fnameR=CPs.getFullFname(database,I,2);
            S=load([fnameL '.mat']);
            CPs_db.check_fun_im(S.cps,I,'L',IszRC);
            S=load([fnameR '.mat']);
            CPs_db.check_fun_im(S.cps,I,'R',IszRC);

        end
    end
    function out= check_xyz(database,startAt)
        if ~exist('startAt','var') || isempty(startAt)
            startAt=1;
        end
        db=dbInfo(database);
        ActrRC=CPs.getAllRC(db.IszRC);
        dire=CPs.getFullDire(database);
        if ~Dir.exist(dire)
            mkdir(dire);
        end
        N=length(db.allImages);
        p=Pr(N,1,'Generating CPsFull');
        IszRC=db.IszRC;

        for i = 1:N
            p.u();
            if i < startAt
                continue
            end
            I=db.allImages(i);
            cpLookup=CPs.loadXYZ(database,I,0);
            str=CPs_db.check_fun_lookup(cpLookup{1}{1},'L','L',IszRC,I);
            p.append_msg(str);
            str=CPs_db.check_fun_lookup(cpLookup{1}{2},'L','R',IszRC,I);
            p.append_msg(str);
            str=CPs_db.check_fun_lookup(cpLookup{2}{1},'R','L',IszRC,I);
            p.append_msg(str);
            str=CPs_db.check_fun_lookup(cpLookup{2}{2},'R','R',IszRC,I);
            p.append_msg(str);
        end
        p.c();
    end

    function out= check_lookup(database,startAt)
        if ~exist('startAt','var') || isempty(startAt)
            startAt=1;
        end
        db=dbInfo(database);
        ActrRC=CPs.getAllRC(db.IszRC);
        dire=CPs.getFullDire(database);
        if ~Dir.exist(dire)
            mkdir(dire);
        end
        N=length(db.allImages);
        p=Pr(N,1,'Generating CPsFull');
        IszRC=db.IszRC;
        for i = 1:N
            p.u();
            if i < startAt
                continue
            end
            I=db.allImages(i);
            cpLookup=CPs.loadLookup(database,I,0,0);
            str=CPs_db.check_fun_lookup(cpLookup{1}{1},'L','L',IszRC,I);
            p.append_msg(str);
            str=CPs_db.check_fun_lookup(cpLookup{1}{2},'L','R',IszRC,I);
            p.append_msg(str);
            str=CPs_db.check_fun_lookup(cpLookup{2}{1},'R','L',IszRC,I);
            p.append_msg(str);
            str=CPs_db.check_fun_lookup(cpLookup{2}{2},'R','R',IszRC,I);
            p.append_msg(str);
            waitforbuttonpress
        end
        p.c();
    end
    function str=check_fun_lookup(itpRC,Anchor,LorR,IszRC,I)
        str='';
        naninds=isnan(itpRC);
        gdinds=mod(itpRC,1)==0.5;
        ng=naninds & gdinds;
        YGd=sum(gdinds(:,1) | naninds(:,1),1);
        bRev=sum(gdinds(:,2) | naninds(:,2),1);
        mmY=Num.minMax(itpRC(:,1));
        mmX=Num.minMax(itpRC(:,2));
        %itpRC(100:105,:)
        if YGd < bRev
            str=[num2str(I) ' ' Anchor ' ' LorR];;
        end
    end
    function check_fun_im(cps,I,LorR,IszRC);
        bSuccess=false(4,1);
        bSuccess(1)=CPs_db.check_fun_fun(cps.AitpRC,'AitpRC',IszRC,I,LorR);
        bSuccess(2)=CPs_db.check_fun_fun(cps.BitpRC,'BitpRC',IszRC,I,LorR);
        bSuccess(3)=CPs_db.check_fun_fun(cps.LitpRC,'LitpRC',IszRC,I,LorR);
        bSuccess(4)=CPs_db.check_fun_fun(cps.RitpRC,'RitpRC',IszRC,I,LorR);
        if any(~bSuccess)
            waitforbuttonpress;
        end
        function bSuccess=check_fun_fun(itpRC,name,IszRC,I,LorR)
            s=sum(abs(mod(itpRC,1))>0==0,1);
            mode=numel(itpRC(:,1))==s;
            bSuccess=false;
            try
                assert(isequal(mode,[1,0]));
                assert(isequal(Num.minMax(itpRC(:,1)),[1 IszRC(1)]));
                bSuccess=true;
            catch
                disp([num2str(I) LorR]);
                disp([ '    ' name]);
                disp([ '      mod : ' num2str(mode)]);
                disp([ '      mm 1: ' num2str(Num.minMax(itpRC(:,1))) ]);
                disp([ '      mm 2: ' num2str(Num.minMax(itpRC(:,2))) ]);
            end

        end
    end
end
methods(Static,Access=private)
end
end
