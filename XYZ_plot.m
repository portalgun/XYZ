classdef XYZ_plot < handle
methods
    function imagesc_xyz_compare(obj)
        RC=[3,4];
        for i =1:3
        for k = 1:2
            rc=[i k+(1*k-1)];
            obj.imagesc_fun('xyz',k,i,rc,RC);
            c=caxis;
            rc=[i (k+1)+(1*k-1)];
            obj.imagesc_fun('itpXYZ',k,i,rc,RC);
            caxis(c);
        end
        end
    end
    function imagesc_itpXYZ(obj)
        Fig.new()

        RC=[3,2];
        c=0;
        for i =1:3
        for k = 1:2
            rc=[i k];
            obj.imagesc_fun('itpXYZ',k,i,rc,RC);
        end
        end
    end
    function imagesc_xyz(obj)
        Fig.new()

        RC=[3,2];
        c=0;
        for i =1:3
        for k = 1:2
            rc=[i k];
            obj.imagesc_fun('xyz',k,i,rc,RC);
        end
        end

    end
    function imagesc_fun(obj,fld,k,xyzNum,rc,RC)
        subPlot(RC,rc(1),rc(2));
        imagesc(obj.(fld){k}(:,:,xyzNum));
        Fig.format('','',[fld obj.LANDR(k)]);
        Fig.formatIm;
        colorbar;
    end

    function scatter_itpXYZ(obj)
        Fig.new()
        plot3(obj.itpXYZ{1}(:,1), obj.itpXYZ{1}(:,2),obj.itpXYZ{1}(:,3),'k.'); hold on;
        plot3(obj.itpXYZ{2}(:,1), obj.itpXYZ{2}(:,2),obj.itpXYZ{1}(:,3),'r.'); hold off;
        ylabel('Y'); xlabel('X'); zlabel('Z');
    end
    function scatter_xyz(obj)
        Fig.new()
        plot3(obj.xyz{1}(:,:,1), obj.xyz{1}(:,:,2),obj.xyz{1}(:,:,3),'b.'); hold on;
        plot3(obj.xyz{2}(:,:,1), obj.xyz{2}(:,:,2),obj.xyz{2}(:,:,3),'m.'); hold off;
        ylabel('Y'); xlabel('X'); zlabel('Z');
    end
end
end
