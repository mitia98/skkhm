 
function map=SKKHM_fun(x,N,sigma,beta,Q)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This code implements the paper:
%%%  Li Q., Mitianoudis N., Stathaki T., "Spatial Kernel K-Harmonic Means 
%%%  Clustering for Multi-spectral Image Thresholding ", 
%%%  IET proceedings on Image Processing, Vol. 1, No. 2, pp. 156-167, 
%%%  June 2007.
%%%
%%% INPUTS:  
%%%           x             Multidimensional input image arranged in a 3D array.
%%%                         x(:,:,i) represents the i-th image    
%%%           N             The number of desired clusters
%%%         sigma           Defines the sigma of the radial basis function (RBF)
%%%                         for the kernel.
%%%         beta            The parameter b controls the magnitude of influence 
%%%                         (when b->0 the influence of neighbourhood is vanishing)
%%%           Q             Denotes a QxQ window neighbourhood around each pixel
%%%                         for the spatial influence.
%%%  For recommended values on sigma, beta, Q, please check the relevant section
%%%  in the paper.
%%% OUTPUTS:
%%%         map             The extracted segmentation map.
%%%
%%%%%%%%%%%%%%%%   Imperial College London 2008   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%    Nikolaos Mitianoudis      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 [S1 S2 S3]=size(x);
 
for i=1:S3
 xr(i,:)=reshape(x(:,:,i),1,S1*S2); 
end

 
 c=random('unif',0,1,S3,N);

%  sigma=5;
%  beta=1.2;
%  Q=3;
 L=1;
 

 
  Qc=50000;dQc=Qc;

 u=random('unif',0,1,N,S1*S2);
 u=u./(ones(N,1)*sum(u));
 w=random('unif',0,1,1,S1*S2);
 w=ones(1,S1*S2);
 F=ones(N,S1*S2);
 
 iter=0;
 map=zeros(S1,S2);

 while dQc>0.01
     iter=iter+1
     
     for i=1:N
         c(:,i)=mean((ones(S3,1)*(u(i,:).*w)).*xr,2)./mean(u(i,:).*w);
         d(i,:)=exp(-sum((xr-c(:,i)*ones(1,S1*S2)).^2,1)/sigma^2);
         dd(i,:)=max(1-d(i,:),0.0001);
         tmp=filter2(ones(Q,Q),1-reshape(u(i,:),S1,S2),'same');
         tmp=L.*exp(-beta.*tmp);
         F(i,:)=reshape(tmp,1,S1*S2);
     end
    
     
     tmp1=sum((d.*F)./(dd.^2));          
     for i=1:N
         u(i,:)=(d(i,:).*F(i,:))./(dd(i,:).^2.*tmp1);
     end
     w=tmp1./(sum(F./dd).^2);
    
  
    
     
     u_max=max(u);du=u-ones(N,1)*u_max;

     for i=1:N
         map(find(du(i,:)==0))=i-1;
     end
    imagesc(reshape(map,S1,S2));
     pause(0.1);



     Qc_n=sum(2*2./sum(F./dd));
     dQc=abs(Qc_n-Qc);
     Qc=Qc_n;
    % QQQ=[QQQ Qc];
     
 end

 map=reshape(map,S1,S2);
