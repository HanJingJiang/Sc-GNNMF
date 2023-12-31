clc,clear all
tic%时间开始
load('LD1.mat') % 查找邻接矩阵五折后的训练集和测试集
load('S_L_D.mat')%查找LncRNA和Disease 的相似性局矩阵
load("Data.mat")

oldLD = LD;
adLD = oldLD >0;

Y=LD_new; %Y=train{1}; train{2}; train{3}; train{4}
%（最佳参数选取：k=90，lmadp=lmadl=0.1;beta1=beta2=0.01;）


%params.k=90;
params.iterate=3000; 
  

params.k = 20;
k=params.k;


values=[0.1];
for i = 1:length(values)
    params.beta1=0.01;
    params.beta2=0.01;
    params.lmadp=values(i);
    params.lmadl=values(i); 
    iterate=params.iterate; 
    lmadl=params.lmadl;  
    lmadd=params.lmadp; 
    beta1=params.beta1;
    beta2=params.beta2; 
    fprintf('k=%d  maxiter=%d  lmadp=%d lmadl=%d  beta1=%d beta2=%d\n', k, iterate,lmadl,lmadd,beta1,beta2);

    [n,m]=size(Y);
    %U=rand(k,n);
    %V=rand(k,m);
    [W,H]=NNDSVD(double(Y),k,2);
    U=W';
    V=H;


  
%save WH_value;
    D_L=zeros(n);
    D_D=zeros(m);
    for i=1:n
        D_L(i,i)=sum(S_L(i,1:n));%D_d(i,i)取值为S_d的第i行之和，非主对角元取值为0;D_c同
    end
    for i=1:m
        D_D(i,i)=sum(S_D(i,1:m));
    end
    L_L=D_L-S_L;
    L_D=D_D-S_D;
    

%——————————————迭代法则————————————————
    fid = fopen( 'RunResult.txt','wt+');
    for step=1:iterate
        U1=U.*((V*Y'+lmadl*U*S_L)./(V*V'*U+lmadl*U*D_L+beta1*U));%迭代生成的U矩阵
       
        V1=V.*((U*Y+lmadd*V*S_D)./(U*U'*V+lmadd*V*D_D+beta2*V));%迭代生成的H矩阵
      
     
        ALA = sum(diag((U*L_L)*U'));%ALA = sum(diag((U*L_L)*U'))=Tr((U*L_L)*U')=sum(sum((U*L_L).*U));;
        BLB = sum(diag((V*L_D)*V'));%BLB = sum(diag((V*L_D)*V'))=Tr((V*L_D)*V')=sum(sum((H*L_D).*V));
        obj = sum(sum((Y-U'*V).^2))+beta1*(sum(sum(U.^2)) )+beta2*(sum(sum(V.^2)))+lmadl*ALA+lmadd*BLB;
        error=max([max(sum((U1-U).^2)),max(sum((V1-V).^2))]);%用列F-范数来保证收敛
      
        fprintf(fid,'%s\n',[sprintf('step = \t'),int2str(step),...
                sprintf('\t obj = \t'),num2str(obj),...
		        sprintf('\t error = \t'), num2str(error)]);
                fprintf('step=%d  obj=%d  error=%d\n',step, obj, error); 
        if error< 10^(-4)
                fprintf('step=%d\n',step);
                break;
        end
        U=U1; V=V1;    %将U1,V1赋值给U,V  
    end
    fclose(fid);
    Y5=U'*V; %5个训练集分别得到5个预测结果

    file_name=sprintf('Y5.csv', lmadl);
    writematrix(Y5,file_name)
    writematrix(U,'U.csv')
    writematrix(V,'V.csv')
end


