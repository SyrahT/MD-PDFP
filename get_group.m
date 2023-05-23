function [B,vd] = get_group(d,type)
% Matrix corresponding to the Fused Lasso
if(type==1)
    vd = d-1;
    bx = zeros([2*(d-1),1]);
    by = zeros([2*(d-1),1]);
    bz = zeros([2*(d-1),1]);
    for i=1:d-1
        idx = 2*i-1;
        bx(idx) = i;
        by(idx) = i;
        bz(idx) = 1;
        idx = 2*i;
        bx(idx) = i;
        by(idx) = i+1;
        bz(idx) = 1;
    end
    B = sparse(bx,by,bz,d-1,d);
elseif(type==2)
    vd = 2*d-1;
    bx = zeros([d+2*(d-1),1]);
    by = zeros([d+2*(d-1),1]);
    bz = zeros([d+2*(d-1),1]);
    for i=1:d-1
        idx = 2*i-1;
        bx(idx) = i;
        by(idx) = i;
        bz(idx) = 1;
        idx = 2*i;
        bx(idx) = i;
        by(idx) = i+1;
        bz(idx) = 1;
    end
    for i=1:d
        bx(2*(d-1)+i) = i+d-1;
        by(2*(d-1)+i) = i;
        bz(2*(d-1)+i) = 1;
    end
    B = sparse(bx,by,bz,2*d-1,d);
elseif(type==3)
    vd = d-1;
    bx = zeros([2*(d-1),1]);
    by = zeros([2*(d-1),1]);
    bz = zeros([2*(d-1),1]);
    for i=1:d-1
        idx = 2*i-1;
        bx(idx) = i;
        by(idx) = i;
        bz(idx) = 1;
        idx = 2*i;
        bx(idx) = i;
        by(idx) = i+1;
        bz(idx) = -1;
    end
    B = sparse(bx,by,bz,d-1,d);
elseif(type==4)
    p1 = 2;
    p2 = 5;
    vd = 2*d-p1-p2;
    bx = zeros([2*(2*d-p1-p2),1]);
    by = zeros([2*(2*d-p1-p2),1]);
    bz = zeros([2*(2*d-p1-p2),1]);
    idx = 0;
    for i=1:p1
        for j=1:p2-1
            idx = idx + 1;
            edge = idx*2-1;
            node = (i-1)*p2 + j;
            bx(edge) = idx;
            by(edge) = node;
            bz(edge) = 1;
            edge = idx *2;
            bx(edge) = idx;
            by(edge) = node + 1;
            bz(edge) = -1;
        end
    end
    for i=1:p1-1
        for j=1:p2
            idx = idx + 1;
            edge = idx*2-1;
            node = (i-1)*p2 + j;
            bx(edge) = idx;
            by(edge) = node;
            bz(edge) = 1;
            edge = idx *2;
            bx(edge) = idx;
            by(edge) = node + p2;
            bz(edge) = -1;
        end
    end
    B = sparse(bx,by,bz,vd,d);
end
end

