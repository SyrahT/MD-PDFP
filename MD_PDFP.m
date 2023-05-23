function [re] = MD_PDFP(x_star, gamma, lam, round,A,y)
% MD-PDFP algorithm
    [p,d] = size(A,[2,3]);
    global W;
    W = get_W(d,4);
    rho = norm(W - ones(d)/d,2);
    
    mu = 1;
    [B,vd] = get_group(p,3);
    
    x0 = zeros([d,p]);
    v0 = zeros([d,vd]);
    
    nablafx0 = zeros([d,p]);
    atb = zeros([d,p]);
    for i=1:d
        Ai = A(:,:,i);
        yi = y(:,i);
        nablafx0(i,:) = Ai'*Ai*x0(i,:)';
        atb(i,:) = Ai'*yi;
    end
    nablafx0  = nablafx0 - atb;
    gt = nablafx0;
    
    K0 = 25;
    ibbt = sparse(1:vd,1:vd,1,vd,vd)- lam*(B*B');
    
    re = [];
    for k=1:round
        energy = mu*norm(B*x0(1,:)',1);
        for i=1:d
            Ai = A(:,:,i);
            yi = y(:,i);
            energy = energy + 0.5*norm(Ai*x0(1,:)'-yi)^2/d;
        end
        re = [re,norm(x0-ones([d,1])*x_star','fro')];
        
        K = K0;
        vt = (x0-gamma*gt)*B'+v0*ibbt;
        
        thresh = mu*gamma/lam;
        vt(vt>thresh) = thresh;
        vt(vt<-thresh) = -thresh;
        vt = fastmix(vt,K,rho);
        xt = x0 - gamma*gt - lam*(vt*B);
        xt = fastmix(xt,K,rho);
        x0 = xt;
        
        nablafxt = zeros([d,p]);
        for i=1:d
            Ai = A(:,:,i);
            nablafxt(i,:) = Ai'*(Ai*x0(i,:)');
        end
        nablafxt = nablafxt - atb;
        gt = gt + nablafxt - nablafx0;
        gt = fastmix(gt,K,rho);
        nablafx0 = nablafxt;
        v0 = vt;
    end
end