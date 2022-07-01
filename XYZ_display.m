classdef XYZ_display  < handle
methods
    function [Lmap,Rmap]=xyz_to_ptch_msk(obj,pointsXYZ,dispORdb)
        [PPxyL,PPxyR]=obj.xyz_to_CPs(pointsXYZ,dispORdb);
        [Lmap,Rmap]=obj.CPs_to_patch_msk(PPxyL,PPxyR,PszXY);
    end
    function [Lmap,Rmap]=CPs_to_patch_msk(obj,PPxyL,PPxyR,PszXY)
        Lmap=Msk.rec(obj.IszRC,PPxyL,PszXY(2),PszXY(1),0);
        Rmap=Msk.rec(obj.db.IszRC,PPxyR,PszXY(2),PszXY(1),0);
    end
    function [PPxyL,PPxyR]=xyz_to_CPs(obj,pointsXYZ,dispORdb)
    % for tracing patch/sample window
        IppXm=obj.(dispORdb).IppXm;
        IppYm=obj.(dispORdb).IppYm;
        IppZm=obj.(dispORdb).IppZm;
        if isempty(obj.(dispORdb).X)
            obj.(dispORdb).getXY();
        end
        obj.X=obj.(dispORdb).X;
        obj.Y=obj.(dispORdb).Y;
        [PPxyL,PPxyR]=back_project(obj,pointsXYZ,IppXm,IppYm,IppZm,X,Y);
    end
%%
    function obj=get_VDisp_xyz(obj,displ,PctrDegXY,PszDegXY,IPDm)
        if isempty(obj.IPDm)
            IPDm=.065;
        end
        pixPerMxy=disp.pixPerMxy;
        prsntSzRC=fliplr(PszDegXY.*disp.degPerMXY);

        prsntCtrXYZm=sign(PctrDegXY).*displ.scrnZm.*tand(abs(PctrDegXY));
        prsntCtrXYZm(3)=displ.scrnZm;

        obj.dispXYZ=obj.get_dispXYZ(ind,prsntCtrXYZm,prsntSzRC,pixPerMxy,IPDm);
    end

    function pointsXYZ=get_dispXYZ(obj,k,prsntCtrXYZm,prsntSzRC,pixPerMxy,IPDm)
        pointsXYZ=reshape(obj.xyz{k},obj.db.nPix,3);
        if isempty(obj.db.X)
            obj.db.getXY();
        end

        [PPxyL,PPxyR]=obj.back_project(pointsXYZ,obj.db.IppXm,obj.db.IppYm,obj.db.X,obj.db.Y);
        obj.disp=disp(prsntCtrXYZm,prsntSzRC,pixPerMxy,IPDm);
        obj.disp.getXY();
        pointsXYZ=XYZ.forward_project(obj,PPxyL,PPxyR,obj.disp.IppXm,obj.disp.IppYm,obj.disp.IppZm,obj.disp.X,obj.disp.Y);
    end

    function [IppXm,IppYm,X,Y]=get_dsp_proj_plane(obj,IppZm,PctrXYZ,PszXY,pixPerM)
        % XXX IN DISP
    end
end
end
