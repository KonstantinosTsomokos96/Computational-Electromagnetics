%Dimensions nx,ny either 100x50 or 50x25
nx=100; %x direction subdivision mesh
ny=50; %y  direction subdivision mesh
N=(nx-1)*(ny-1); %number of unknowns
h=1/ny/100; %step size(square grid subdivision)
%Population of sparse 5-diagonal matrix
A=4*diag(ones(N,1))-diag(ones(N-1,1),1)-diag(ones(N-1,1),-1)-diag(ones(N-(nx-1),1),-(nx-1))-diag(ones(N-(nx-1),1),nx-1);
k=nx-1;
for i=1:(ny-2)
    A(k,k+1)=0; 
    A(k+1,k)=0; 
    k=k+(nx-1); 
end
%TE boundary nodes(Neumann condition)
m=1;
for i=1:2 %first and last line of mesh
A(m,m)=2; %1st corner node of line
  for m=m+1:m+(nx-3) 
     A(m,m)=3;
  end
A(m+1,m+1)=2;%2nd corner node of line
m=1+(nx-1)*(ny-2);
end
m=nx; %near edge nodes
for i=1:(ny-3)
    A(m,m)=3;
    A(m+nx-2,m+nx-2)=3;
    m=m+nx-1;
end
%wavenumber k cutoff - theoretical analysis
k_th = [ 0 157.08 314.15 314.15 351.23]' ; 
%wavelength
lambda = eigs(A,5,'smallestabs') ;
%k cutoff - numerical method
k_num = sqrt(lambda)/h ;
%errors
errors = abs(k_th - k_num)./k_th *100;
fprintf(' TE_10 error: %f%% \n TE_01 error: %f%% \n TE_20 error: %f%% \n TE_11 error: %f%% \n',errors(2),errors(3),errors(4),errors(5))
%eigenvalues and eigenvectors of sparse matrix
[V,D]=eig(A); 
[d,ind] = sort(diag(D)); 
V = V(:,ind);
%Matrix of each mode , Magnetic field
%TE10
H1=zeros(ny+1,nx+1);
H1(2:ny,2:nx)=reshape(V(:,2),[nx-1,ny-1])' ;
H1(2:ny,1)=H1(2:ny,2);
H1(2:ny,nx+1)=H1(2:ny,nx);
H1(1,1:nx+1)=H1(2,1:nx+1);
H1(ny+1,1:nx+1)=H1(ny,1:nx+1);
%TE20
H2=zeros(ny+1,nx+1);
H2(2:ny,2:nx)=reshape(V(:,3),[nx-1,ny-1])' ;
H2(2:ny,1)=H2(2:ny,2);
H2(2:ny,nx+1)=H2(2:ny,nx);
H2(1,1:nx+1)=H2(2,1:nx+1);
H2(ny+1,1:nx+1)=H2(ny,1:nx+1);
%TE01
H3=zeros(ny+1,nx+1);
H3(2:ny,2:nx)=reshape(V(:,4),[nx-1,ny-1])' ;
H3(2:ny,1)=H3(2:ny,2);
H3(2:ny,nx+1)=H3(2:ny,nx);
H3(1,1:nx+1)=H3(2,1:nx+1);
H3(ny+1,1:nx+1)=H3(ny,1:nx+1);
%TE11
H4=zeros(ny+1,nx+1);
H4(2:ny,2:nx)=reshape(V(:,5),[nx-1,ny-1])' ;
H4(2:ny,1)=H4(2:ny,2);
H4(2:ny,nx+1)=H4(2:ny,nx);
H4(1,1:nx+1)=H4(2,1:nx+1);
H4(ny+1,1:nx+1)=H4(ny,1:nx+1);
%contour plots
x=0:h*100:2; 
y=0:h*100:1;
[X,Y]=meshgrid(x,y);
subplot(2,2,1)
contour(X,Y,H1,30),title('Magnetic Field , TE10 mode'),xlabel('x'),ylabel('y')
colorbar
subplot(2,2,2)
contour(X,Y,H2,30),title('Magnetic Field , TE20 mode'),xlabel('x'),ylabel('y')
colorbar
subplot(2,2,3)
contour(X,Y,H3,30),title('Magnetic Field , TE01 mode'),xlabel('x'),ylabel('y')
colorbar
subplot(2,2,4)
contour(X,Y,H4,30),title('Magnetic Field , TE11 mode'),xlabel('x'),ylabel('y')
colorbar
