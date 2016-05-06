%Initialization
Y=X4I9DF0timeseries;
Sigma0=corr(Y');
U=corr(Y');
V=corr(Y');
A=eye(167);
B=eye(167);

N=1200;
pred_mean=cell(1,N);
pred_cov=cell(1,N);
filt_mean=cell(1,N+1);
filt_cov=cell(1,N+1);
smooth_mean=cell(1,N+1);
smooth_cov=cell(1,N+1);

%EM
for rep=1:100
%Estep
filt_mean{1}=zeros(167,1);
filt_cov{1}=Sigma0;
%%Kalman filter
for i=1:N
   pred_mean{i}=A*filt_mean{i};
   pred_cov{i}=A*filt_cov{i}*A'+U;

   filt_cov{i+1}=inv(inv(pred_cov{i})+B'*(V\B));
   filt_mean{i+1}=pred_mean{i}+filt_cov{i+1}*B'*(V\(Y(:,i)-B*pred_mean{i}));
end


%smoothing
L=cell(1,N);
for i=1:N
    L{i}=filt_cov{i}*A'/(pred_cov{i});
end

smooth_mean{N+1}=filt_mean{N+1};
smooth_cov{N+1}=filt_cov{N+1};
for i=N:-1:1
   smooth_mean{i}=filt_mean{i}+L{i}*(smooth_mean{i+1}-pred_mean{i});
   smooth_cov{i}=filt_cov{i}+L{i}*(smooth_cov{i+1}-pred_cov{i})*L{i}';
end
  
xxt=cell(1,N);
xx_1t=cell(1,N);

for i=1:N
    xxt{i}=smooth_cov{i+1}+smooth_mean{i+1}*smooth_mean{i+1}';
    temp=L{i}*smooth_cov{i+1}+(filt_mean{i}+L{i}*(smooth_mean{i+1}-pred_mean{i})-smooth_mean{i})*smooth_mean{i+1}';
    xx_1t{i}=temp+smooth_mean{i+1}*smooth_mean{i}';
end

%Mstep
Sigma0=smooth_cov{1}+smooth_mean{1}*smooth_mean{1}';
 
temp1=zeros(167,167);
temp2=Sigma0;
for i=1:N
    temp1=temp1+xx_1t{i};
end
for i=1:(N-1)
    temp2=temp2+xxt{i};
end
A=temp1/(temp2); 


temp1=zeros(167,167);
temp2=zeros(167,167);
for i=1:N
    temp1=temp1+Y(:,i)*smooth_mean{i+1}';
    temp2=temp2+xxt{i};
end
B=temp1/(temp2);


U=xxt{1}-A*xx_1t{1}'-xx_1t{1}*A'+A*Sigma0*A';
for i=2:N
    U=U+xxt{i}-A*xx_1t{i}'-xx_1t{i}*A'+A*xxt{i-1}*A';
end
U=U/N;
U=(U+U')/2;

V=zeros(167,167);
for i=1:N
    V=V+Y(:,i)*Y(:,i)'-B*smooth_mean{i+1}*Y(:,i)'-Y(:,i)*smooth_mean{i+1}'*B'+B*xxt{i}*B';
end
V=V/N;
V=(V+V')/2;
end

