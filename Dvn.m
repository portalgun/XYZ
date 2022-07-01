classdef Dvn < handle
methods(Static)
    function Ximg=nInImg(DVN,k,cpLookup,IszRC,kernSz,thresh);
        [k,nk]=CPs.getK(k);

        DVN{1}=uint16(DVN{1});
        DVN{2}=uint16(DVN{2});

        [ROWS,COLS]=Range.getConvRC(IszRC,kernSz);
        AitpRC0=CPs.getAitpRC0(kernSz);
        Ximg=nan(IszRC);

        for i = ROWS % ROWS
        for j = COLS%  COLS
            ActrRC=[i,j];
            [AitpInd,BitpInd]=CPs.ActrRCToABRecInd(ActrRC,AitpRC0,k,cpLookup,IszRC);

            if any(CPs.getIndOutOfRangeInd([AitpInd; BitpInd],IszRC))
                continue
            end

            Amap=reshape(DVN{k}(AitpInd),kernSz)';
            try
                Bmap=reshape(DVN{nk}(BitpInd),kernSz)';
            catch
                continue
            end

            %% UNIQUE
            Acnt=max(Msk.nRunLengthFast(Amap,2,[],thresh));
            Bcnt=max(Msk.nRunLengthFast(Bmap,2,[],thresh));
            Ximg(i,j)=Acnt+Bcnt;
            %%

        %    Ximg(i,j)=Acnt+Bcnt;
        end
        end
    end
    function Ximg=nRegInImg(DVN,k,cpLookup,IszRC,kernSz,thresh);
        [k,nk]=CPs.getK(k);

        DVN{1}=uint16(DVN{1});
        DVN{2}=uint16(DVN{2});

        [ROWS,COLS]=Range.getConvRC(IszRC,kernSz);
        AitpRC0=CPs.getAitpRC0(kernSz);
        Ximg=nan(IszRC);

        for i = ROWS % ROWS
        for j = COLS%  COLS
            ActrRC=[i,j];
            [AitpInd,BitpInd]=CPs.ActrRCToABRecInd(ActrRC,AitpRC0,k,cpLookup,IszRC);

            if any(CPs.getIndOutOfRangeInd([AitpInd; BitpInd],IszRC))
                continue
            end

            Amap=reshape(DVN{k}(AitpInd),kernSz)';
            try
                Bmap=reshape(DVN{nk}(BitpInd),kernSz)';
            catch
                continue
            end


            %% UNIQUE
            ACC = bwconncomp(Amap);
            BCC = bwconncomp(Bmap);
            Ximg(i,j)=ACC.NumObjects + BCC.NumObjects;
            %%
        end
        end
    end
    function Ximg=magInImg(DVN,k,cpLookup,IszRC,kernSz,thresh);
        [k,nk]=CPs.getK(k);

        DVN{1}=uint16(DVN{1});
        DVN{2}=uint16(DVN{2});

        [ROWS,COLS]=Range.getConvRC(IszRC,kernSz);
        AitpRC0=CPs.getAitpRC0(kernSz);
        Ximg=nan(IszRC);

        for i = ROWS % ROWS
        for j = COLS%  COLS
            ActrRC=[i,j];
            [AitpInd,BitpInd]=CPs.ActrRCToABRecInd(ActrRC,AitpRC0,k,cpLookup,IszRC);

            if any(CPs.getIndOutOfRangeInd([AitpInd; BitpInd],IszRC))
                continue
            end

            Amap=reshape(DVN{k}(AitpInd),kernSz)';
            try
                Bmap=reshape(DVN{nk}(BitpInd),kernSz)';
            catch
                continue
            end


            %% UNIQUE
            Ximg(i,j)=sum(Amap,'all')+sum(Bmap,'all');
            %%

        end
        end
    end
end
end
