%Dimensions nx,ny either 100x50 or 50x25
nx=100; %x direction subdivision mesh
ny=50; %y  direction subdivision mesh
h=1/ny/100;%step size(square grid subdivision)
N=(nx-1)*(ny-1); %number of unknowns
%Population of sparse 5-diagonal matrix
A=4*diag(ones(N,1))-diag(ones(N-1,1),1)-diag(ones(N-1,1),-1)-diag(ones(N-(nx-1),1),-(nx-1))-diag(ones(N-(nx-1),1),nx-1);
k=nx-1;
for i=1:(ny-2) 
    A(k,k+1)=0; 
    A(k+1,k)=0; 
    k=k+(nx-1); 
end
%wavenumber k cutoff - theoretical analysis
k_th = [ 351.23 444.27 566.34 647.63 ]' ;
%wavelength
lambda = eigs(A,4,'smallestabs');
%wavenumber k cutoff - numerical method
k_num = sqrt(lambda)/h ;
%errors
errors = abs(k_th - k_num)/k_th *100 ;
fprintf(' TM_11 error: %f%% \n TM_21 error: %f%% \n TM_31 error: %f%% \n TM_12 error: %f%%\n\n',errors(1),errors(2),errors(3),errors(4))
%eigenvalues and eigenvectors of sparse matrix
[V,D]=eig(A); 
[d,ind] = sort(diag(D)); %sorted eigenvalues in ascending order 
V = V(:,ind) ;
%Matric of each mode , Electric field
%TM11
E1=zeros(ny+1,nx+1);
E1(2:ny,2:nx)=reshape(V(:,1),nx-1,ny-1)' ;
%TM21
E2=zeros(ny+1,nx+1);
E2(2:ny,2:nx)=reshape(V(:,2),nx-1,ny-1)' ;
%TM31
E3=zeros(ny+1,nx+1);
E3(2:ny,2:nx)=reshape(V(:,3),nx-1,ny-1)' ;
%тл12
E4=zeros(ny+1,nx+1);
E4(2:ny,2:nx)=reshape(V(:,4),nx-1,ny-1)' ;
%contour plots
x=0:h*100:2;
y=0:h*100:1;
[X,Y]=meshgrid(x,y);
%TM11
subplot(2,2,1)
contour(X,Y,E1,20),title('Electric Field, TM11 mode'),xlabel('x'),ylabel('y')
colorbar
%TM21
subplot(2,2,2)
contour(X,Y,E2,20),title('Electric Field, TM21 mode'),xlabel('x'),ylabel('y')
colorbar
%TM31
subplot(2,2,3)
contour(X,Y,E3,20),title('Electric Field, TM31 mode'),xlabel('x'),ylabel('y')
colorbar
%TM12
subplot(2,2,4)
contour(X,Y,E4,20),title('Electric Field, TM12 mode'),xlabel('x'),ylabel('y')
colorbar

