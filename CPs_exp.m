classdef CPs_exp < handle
methods(Static)
    function cpRC=lookupRC2(AitpRC,k,cpLookup,IszRC)
        k=CPs.getK(k);
        N=size(cpLookup{k},1)-1;
        n=size(AitpRC,1);
        IszC=cpLookup{k}(end,1);
        medD=cpLookup{k}(end,2);

        inds = AitpRC(:,2) + (AitpRC(:,1)-1)*IszC;
        cpRC=cpLookup{k}( inds ,:);
    end
    function cpRC=lookupRC3(AitpRC,k,cpLookup,IszRC)
        k=CPs.getK(k);
        N=size(cpLookup{k},1)-1;
        n=size(AitpRC,1);
        IszC=cpLookup{k}(end,1);
        medD=cpLookup{k}(end,2);

        inds = AitpRC(:,2)-medD + ((AitpRC(:,1))-1)*IszC;
        cpRC=nan(n,2);
        gdInds=inds <= N & inds > 0;
        cpRC(gdInds,:)=cpLookup{k}( inds(gdInds) ,:);
    end
    function [cpLookup,exitflag1,exitflag2]=toLookup(CPs_L_AB,CPs_R_AB,db,ActrRC)
        cpLookup=cell(2,1);
        [cpLookup{1},exitflag1]=CPs.toLookupLorR(CPs_L_AB,'L',db,ActrRC);
        [cpLookup{2},exitflag2]=CPs.toLookupLorR(CPs_R_AB,'R',db,ActrRC);
    end

    function [cpLookup,exitflag]=toLookupLorR(CPs_A_AB,LorR,db,ActrRC)
        exitflag=false;
        AitpRC=CPs_A_AB{1};
        AitpRCO=CPs_A_AB{1};
        BitpRC=CPs_A_AB{2};
        A=ActrRC;
        B=BitpRC;
        indNan2=isnan(AitpRC(:,2));
        indNan1=isnan(AitpRC(:,1));

        %%% FILL
        % FILL IN FIRST ROW, DOES NOT EVER CHANGE
        AitpRC(indNan1,1)=ActrRC(indNan1,1);

        %%% SORT
        [A,B]=CPs_exp.sortBToNearestPlaceByA(AitpRC,BitpRC,ActrRC,db.IszRC);

        %%% CHECK
        if nansum(round(A)-ActrRC,'all')>0.1
            exitflag=true;
            cpLookup=[];
            return
        end
        assert(nansum(round(A)-ActrRC,'all')<0.1);
        if ~isequal(size(B),size(ActrRC))
            size(B)
            size(ActrRC)
            error('sizes do not match');
        end

        if LorR=='L'
            cpLookup{1}=A;
            cpLookup{2}=B;
        elseif LorR=='R'
            cpLookup{1}=B;
            cpLookup{2}=A;
        end

    end
    function [ANew,BNew]=sortBToNearestPlaceByA(A,B,Ideal,IszRC)
        N=prod(IszRC);
        rA=round(A);
        AInds=Sub.toInd(rA,IszRC);
        IdealInds=transpose(1:N);


        D=sqrt(sum(A(:,2)-rA(:,2)).^2);

        % OUT OF RANGE GETS PUT INTO LAST BIN
        badInds=isnan(AInds) | AInds > N | AInds < 1;
        AInds(badInds)=[];
        A(badInds,:)=[];
        B(badInds,:)=[];

        Bins=1:N;
        Edges=[0:N]+0.5;

        % slots where each AInd should go
        binInds=discretize(AInds,Edges);

        %[binInds]=sort(binInds);
        hc=histcounts(binInds,Edges);
        maxc=max(hc);


        ANew=nan(N,2,1);
        BNew=nan(N,2,1);
        %for i = 1:5
        %gdInd=ismember(binInds,Bins(hc==1));
        %sum(hc<1)
        %sum(hc>1)

        method=1;
        if method==1
            gdInd=ismember(binInds,Bins(hc<5));
            gdBins=binInds(gdInd);

            ANew(gdBins,:,:)=A(gdInd,:);
            BNew(gdBins,:,:)=B(gdInd,:);

        elseif method==2
            gdInd=ismember(binInds,Bins(hc==1));
            gdBins=binInds(gdInd);
            ANew(gdBins,:,:)=A(gdInd,:);


            BinPool=Bins(hc>1);
            uBinPool=unique(BinPool);

            % TODO
            %%  SLOWW
            for i=1:length(uBinPool)
                gdInds=find(ismember(binInds,uBinPool(i)));
                gdBins=binInds(gdInds(1));

                d=abs(round(A(gdInds,2))-A(gdInds,2));
                sel=d==min(d);
                Anew(gdBins,:)=A(gdInds(sel),:);
                Bnew(gdBins,:)=B(gdInds(sel),:);
            end

        end


    end
    function cpLooup=toLookupLorR2(CPs_A_AB,LorR,db,ActrRC)
        AitpRC=CPs_A_AB{1};
        AitpRCO=CPs_A_AB{1};
        BitpRC=CPs_A_AB{2};
        A=ActrRC;
        B=BitpRC;
        indNan2=isnan(AitpRC(:,2));
        indNan1=isnan(AitpRC(:,1));
        % FILL IN FIRST ROW, DOES NOT EVER CHANGE
        %%% FILL
        AitpRC(indNan1,1)=ActrRC(indNan1,1);

        %%% SORT
        bSort=false;
        if bSort
            inds=find(~indNan2);
            [~,order]=sortrows(AitpRC(~indNan2,:));
            NanIndNan=inds(order);
            AitpRC(NanIndNan,:)=AitpRC(~indNan2,:);
        end
        %%%

        d=AitpRC(:,2)-ActrRC(:,2);
        r=d-round(d);
        medD=nanmedian(d);

        rndMedD=round(medD);

        inds=find(~indNan2);
        [~,order]=sortrows(AitpRC(~indNan2,:));
        NanIndNan=inds(order);

        if rndMedD > 0
             indNan1=find(isnan(AitpRC(:,1)))+rInd;
             indNan2=find(isnan(AitpRC(:,2)))+rInd;
             indNan1(indNan1 < 1 | indNan1 > size(AitpRC,1))=[];
             indNan2(indNan2 < 1 | indNan2 > size(AitpRC,1))=[];
        end
        A(indNan1,:)=nan;
        A(indNan2,:)=nan;

        nans=nan(abs(rInd),2);
        if rndMedD > 0
            A(1:rInd,:)=nan;

            e=size(A,1);
            inds=e-rInd+1:e;

            B=[nans; B];
            B(inds,:)=[];

            r=[nans(:,1); r];
            r(inds,:)=[];
        elseif rndMedD < 0
            inds=1:rInd;

            A(end-rInd+1:end,:)=nan;
            B=[B; nans];
            B(inds,:)=[];

            r=[nans(:,1); r];
            r(inds,:)=[];
        end
        A(:,2)=A(:,2)+r;
        %assert(nansum(ActrRC(:,2)-(A(:,2)-r))<1);
        %assert(nansum(ActrRC(:,1)-(A(:,1)))<1);


        if ~isequal(size(B),size(ActrRC))
            size(B)
            size(ActrRC)
            error('sizes do not match');
        end
        if LorR=='L'
            cpLookup{1}=A;
            cpLookup{2}=B;
        elseif LorR=='R'
            cpLookup{1}=B;
            cpLookup{2}=A;
        end



        %AitpRC(indNan1,1)=IitpRC(indNan1,1);

        %AitpRC(:,2)=AitpRC(:,2)-d;

        %% SORT IGNORE NANS
        %cpLookup(S.NanIndNan,:)=cpLookup(~S.indNan2,:);


        %cpLookup(end+1,:)=[db.IszRC(2) fmedD];
    end
    function [cpLookup,inds]=toLookupS(CPs_L_AB,CPs_R_AB,IszRC,xdiff,ydiff)
        cpLookup=cell(2,1);
        [cpLookup{1}]=CPs.toLookupLorRS(CPs_L_AB,'L',IszRC,xdiff,ydiff);
        [cpLookup{2}]=CPs.toLookupLorRS(CPs_R_AB,'R',IszRC,xdiff,ydiff);
    end

    function [cpLookup]=toLookupLorRS(CPs_A_AB,LorR,IszRC,xdiff,ydiff)
        % NOT WORKING COMPLELTEY
        AitpRC=CPs_A_AB{1};
        BitpRC=CPs_A_AB{2};


        %% GET X-OFFSET
        [medD,d]=CPs.getXOffset(AitpRC,IszRC);
        fmedD=round(medD);
        dn=d-fmedD;


        cpLookup=BitpRC;
        cpLookup(:,2)=cpLookup(:,2)-dn;

        %IitpRC=CPs.getAllRC(IszRC);

        %% CORRECTIONS
        %SO1 = AitpRC(:,2) > BitpRC(:,2);
        %SO2 = AitpRC(:,2) >= cpLookup(:,2);


        %oInd=SO1 ~= SO2;
        %cpLookup(oInd,2)=nan;


        cpLookup(end+1,:)=[IszRC(2) fmedD];

    end
end
end
