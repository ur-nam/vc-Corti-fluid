function Nd = update3D(Nd, Uc)
% %   function Nd = update3D(Nd, Uc)
% % -----------------------------------------------------------------------
% % 
% %   Program unit: update3D
% %   Description: retrieve nodal displacement from displacement vector
% %
% %   First written: 2009.03.25
% %   Last modified:
% %
% %   by Jong-Hoon Nam
% %   University of Rochester
% %    
    
    NDOF = 6;
    
    tdof = length(Nd.mapF2R);
    rdof = length(Uc);
    
    Uf=zeros(tdof,1);    
    Uf(Nd.mapR2F(1:rdof)) = Uc;
    

    
    dx=Uf(1:NDOF:tdof);    
    dy=Uf(2:NDOF:tdof);
    dz=Uf(3:NDOF:tdof);
    rx=Uf(4:NDOF:tdof);
    ry=Uf(5:NDOF:tdof);
    rz=Uf(6:NDOF:tdof);
  
    nx = size(Nd.X,1);
    ny = size(Nd.X,2);
 
    for jj=1:ny,
        
        kk = (1:nx)' + (jj-1)*nx;
        
        Nd.dx(:,jj) = dx(kk);
        Nd.dy(:,jj) = dy(kk);
        Nd.dz(:,jj) = dz(kk);
        Nd.rx(:,jj) = rx(kk);
        Nd.ry(:,jj) = ry(kk);
        Nd.rz(:,jj) = rz(kk);

    end
    
end