I=imread('Barbara.tif');
%I=rgb2gray(I);
I=I(1:510,1:510);
[r1,c1]=size(I);
m=fix(r1/5);
n=fix(c1/5);
I=gpuArray(I);
 b=5;
%imshow(I);
%读入要嵌入的秘密图像，黑白图像
J1=imread('hei1.jpg');
J1=imresize(J1,[b,10201]);
J2=rgb2gray(J1);
J=im2bw(J2);
[r2,c2]=size(J);
J=im2uint8(J);%加噪ff声，不能是logical类型
J=imnoise(J,'gaussian',0.02);
J=imnoise(J,'poisson');
J=im2bw(J);
% J=gpuArray(J);


jg=zeros(11,3);
...........产生wd个随机序列............

%b=11;
wd=2^b;
rng('shuffle');
p=floor(255*rand([wd,510,510]));


..........先确定每次要嵌入的比特数，然后循环进行嵌入，在此循环里边计算每次的密钥值，然后根据密钥值的不同，选择不同的随机串进行异或......
zc=m*n*b;    %载体图像所需要嵌入的比特数
mc=r2*c2;    %秘密图像拥有的比特数



tic;

    %先把待嵌入的秘密图像依次存入一个一维向量
    JJ=zeros(1,m*n*b);
    %此处转置
    J=J';
    JJ(1,1:r2*c2)=J(:);
    JJ=reshape(JJ,b,[]);
    [j_c,j_r]=size(J);
    x=zeros(b,1);
    for i=1:b
        x(i,1)=2^(i-1);
    end
    x=repmat(x,1,10404);
    y=arrayfun(@mtimes,JJ,x);
    y=num2cell(y,1);
    F=@(x)sum(x{:});
    y_=arrayfun(F,y);
    y_=y_+1;
    %根据每一块对应的y_值生成505*505大小的随机矩阵
    mask=zeros(510,510);
    maxindex=length(y_);
    
    GetMask5(maxindex,p,y_,mask,wd);
%     for index=1:maxindex;
%         取出第一段秘密信息对应的十进制数
%         num_=y_(1,index)+1;
%         q=squeeze(p(num_,:,:));
%         r_index=ceil(index/102)*5-4;
%         c_index_=mod(index-102*fix(index/102),102);
%         if c_index_==0
%             c_index_=102;
%         end
%         c_index=c_index_*5-4;
%         mask(r_index:r_index+4,c_index:c_index+4)=q(r_index:r_index+4,c_index:c_index+4);
%     end
    
%     for i=2:102;
%         for j=2:102;
%             for index=1:maxindex;
%             x1=(i-1)*5-1;
%             x2=i*5;
%            y1=(j-1)*5-1;
%             y2=j*5;
%              num_=y_(1,index)+1;
%         q=squeeze(p(num_,:,:));
%     
%             mask(x1:x2,y1:y2)=q(x1:x2,y1:y2);
%             end
%         end
%     end
%  


% 让开第一行和第一列进行加密
  mmask=mask(6:510,6:510);
  Ip=I(6:510,6:510);

Ijp=arrayfun(@bitxor,Ip,mmask);%加密完成6:510,6:510部分

Ij=I;
Ij(6:510,6:510)=Ijp;
% 

% 对第一行第一列使用选定的秘钥进行加密
key1=p(17,:,:); 
key1=uint8(key1);
key1=reshape(key1,510,510);
for i=1:102;
   i1=(i-1)*5+1;
   i2=i*5;
   Ij(1:5,i1:i2)=bitxor(key1(1:5,i1:i2),Ij(1:5,i1:i2));
end

key2=p(2,:,:); 
key2=uint8(key2);
key2=reshape(key2,510,510);
for i=2:102;
   i1=(i-1)*5+1;
   i2=i*5;
   Ij(i1:i2,1:5)=bitxor(key2(i1:i2,1:5),Ij(i1:i2,1:5));
end


subplot(221);imshow(I);title('原始图像');
subplot(222);imshow(Ij);title('加密图像');


toc;
...............把加密后的图像分别与每一个随机串进行异或，并求出每次的相邻距离之和...................    

tic;

%先用已知key解密第一行第一列
Ijj=Ij;

for i=1:102;
   i1=(i-1)*5+1;
   i2=i*5;
   Ij1r(1:5,i1:i2)=bitxor(key1(1:5,i1:i2),Ij(1:5,i1:i2));
end
Ijj(1:5,:)=Ij1r(1:5,:);


for i=2:102
    i1=(i-1)*5+1;
    i2=i*5;
    Ij1r(i1:i2,1:5)=bitxor(key2(i1:i2,1:5),Ij(i1:i2,1:5));
end
Ijj(:,1:5)=Ij1r(:,1:5);




%1024个随机矩阵与加密图像进行异或

Ij1=permute(repmat(Ij,1,1,wd),[3,1,2]);  %把加密后图像Ij复制1024遍 再排列成1024*510*510
Ij_=arrayfun(@bitxor,Ij1,uint8(p));      %加密图像与秘钥集异或后的1024*510*510图片
Ij__=gather(Ij_);
%Ij__=Ij__(:,6:510,6:510);
Ijjrc=permute(repmat(Ijj,1,1,wd),[3,1,2]);%解密完第一行和第一列后的Ijj复制1024遍

Ijj_=gather(Ijj);
% for i=1:1024;
%     i
%     Ijrep(i,6:510,6:510)=Ij__(i,:,:);
%     Ijrep(i,1:5,1:5)=Ijj_(1:5,1:5);
% end

%Ijj为部分解密的加密图形
%imshow(Ijj_);
Ijrep=permute(repmat(Ijj,1,1,wd),[3,1,2]);%Ijrep为第一行第一列解密完成剩余部分为加密图像的1024个图像
% %把Ijrep的1024个6:510,6:510替换成Ijbit_的6:510,6:510
for i=1:2^b
    i
    Ijrep(i,6:510,6:510)=Ij__(i,6:510,6:510);
end
% 
% Ijbit_=gather(Ijbit);
% Ijj_=gather(Ijj);
% 
Ijrep_=gather(Ijrep);
%此时Ijrep_中有1024张图 每张图第一行第一列为解密完成部分 6:510,6:510为加密图像与秘钥集异或后图像

% Ijrep(1:5,:)=Ijj_(1:5,:);
% imshow(Ijrep);





%------Ijrep_为解出部分(1:5,1:5)与异或后部分(6:510,6:510)拼成的1024*510*510图形

% Ijjrep=repmat(Ijj,1,1,1024);
% Ijjrep=permute(Ijjrep,[3,1,2]);%Ijjrep中存储1024个解好第一行第一列并且剩余部分为加密图形的图片
% Ijjrep_=gather(Ijjrep);


%此时Ij__中存储着加密后图像Ij的（6:510,6:510）与wd个p的（6:510,6:510）异或后的结果

%计算相关性应加入已经加入已经解好的第一行第一列的块

sum1=zeros(wd,102,102);
for index_=1:wd
    
   sum0=zeros(102,102);
   for i=2:102
       
       for j=2:102
           
            x1=5*(i-1)+1;
            x5=5*i;
            y1=5*(j-1)+1;
            y5=5*j;
            %一个5*5的块只能计算4*5=20次
            d1=zeros(30,1);
            count1=1;
            %用第一行的块作参考 提高准确性
            for i1=x1-1:x5
                for j1=y1-1:(y5-1)
                    d1(count1)=double(int16(Ijrep_(index_,i1,j1+1))-int16(Ijrep_(index_,i1,j1)));
                    count1=count1+1;
                end
            end
            
            d2=zeros(30,1);
            count2=1;
            for j2=y1-1:y5
                for i2=x1-1:(x5-1)
                    d2(count2)=double(int16(Ijrep_(index_,i2+1,j2))-int16(Ijrep_(index_,i2,j2)));
                    count2=count2+1;
                end
            end
            
            s=0;
            for i3=1:30
                s=s+abs(d1(i3))+abs(d2(i3));
            end
            sum0(i,j)=s;
        end
    end
   sum1(index_,:,:)=sum0;
 %sum11(index_,:,:)=sum1(index_,2:102,2:102);
end
%sum1中存储着1024个102*102的数据
%去掉sum1中第一行第一列

toc
%     fun=@(block_struct) abs(block_struct.data(1,2)-block_struct.data(1,1))+abs(block_struct.data(1,3)-block_struct.data(1,2))+abs(block_struct.data(2,2)-block_struct.data(2,1))+abs(block_struct.data(2,3)-block_struct.data(2,2))+abs(block_struct.data(3,2)-block_struct.data(3,1))+abs(block_struct.data(3,3)-block_struct.data(3,2))+abs(block_struct.data(2,1)-block_struct.data(1,1))+abs(block_struct.data(2,2)-block_struct.data(1,2))+abs(block_struct.data(2,3)-block_struct.data(1,3))+abs(block_struct.data(3,1)-block_struct.data(2,1))+abs(block_struct.data(3,2)-block_struct.data(2,2))+abs(block_struct.data(3,3)-block_struct.data(2,3));
%     result=blockproc(squeeze(Ij_(i,:,:)),[3,3],fun);
%     sum1(i,:,:)=result;


%分别找出每个块对应的最小和，将随机串与加密图像的对应块进行异或，从而解出正确的原图像

%sum1中存储1024个102*102 的数字


 %[Y,U]=max(sum1)%返回行向量Y和U，Y向量记录A的每列的最大值，U向量记录每列最大值的行号。
% [~,mk]=min(sum1,[],1);
% mk=squeeze(mk);
% mk=mk';
% mk=mk(:);
% mask_=zeros(505,505);
% GetMask5(maxindex,p,mk,mask_,wd);

%  for index__=1:n*m
%         r_index_=ceil(index__/170);
%         c_index_=mod(index__-170*fix(index__/170),170);
%         if c_index_==0
%             c_index_=170;
%         end
%         r_index=r_index_*3-2;
%         c_index=c_index_*3-2;
%         num_=mk(r_index_,c_index_);
%         q=squeeze(p(num_,:,:));
%         mask_(r_index:r_index+2,c_index:c_index+2)=q(r_index:r_index+2,c_index:c_index+2);
%  end

% Ijq=Ij(6:510.6:510);
% mask_q=mask_(6:510.6:510);
%example Ijj=arrayfun(@bitxor,Ij,mask_);


 
% for i=1:101;
%     for j=1:101
%        A=sum1(:,i,j);
%        
%        [B,num]=min(A);
%        fmin(i,j)=num;
%     end
% end

%用101*101个秘钥序号生成一张505*505的mask_


% for index__=1:101*101
%         r_index_=ceil(index__/170);
%         c_index_=mod(index__-170*fix(index__/170),170);
%         if c_index_==0
%             c_index_=170;
%         end
%         r_index=r_index_*3-2;
%         c_index=c_index_*3-2;
%         num_=mk(r_index_,c_index_);
%         q=squeeze(p(num_,:,:));
%         mask_(r_index:r_index+2,c_index:c_index+2)=q(r_index:r_index+2,c_index:c_index+2);
%  end
tic;
%Ijq=Ij(6:510,6:510);%取加密图像的剩余部分（除去第一行 第一列）
%找出sum11中每一块对应的秘钥序号
[qq,mk]=min(sum1,[],1);
mk=reshape(mk,102,102);
%此处转置
mk=mk';
mask_=zeros(510,510);
GetMask5(maxindex,p,mk,mask_,wd);
mmask_=mask_(6:510,6:510);
% mk=mk';
% mk=mk(:);%把mk排成一维的列向量
%A=Ijjrep(6,:,:);

%  
%  
% for i=1:102;
%     i
%     for j=1:102;
%         x1=(i-1)*5+1;
%         x2=i*5;
%         y1=(j-1)*5+1;
%         y2=j*5;
%         num_=mk(i,j);
%         q=p(num_,:,:);
%         q=reshape(q,510,510);
%     %  q=q(6:510,6:510);
%         mask_(x1:x2,y1:y2)=q(x1:x2,y1:y2);
%     end
% end

%加密图像与拼成的mask_异或
Ijjp=arrayfun(@bitxor,Ijp,mmask_);
Ijj(6:510,6:510)=Ijjp;
subplot(223);imshow(Ijj);title('解密图像');
toc

%计算解密后的图像的错误率
% z=abs(int16(Ijj)-int16(Ipp));
[r3,c3]=size(I);
pc=r3*c3;
ch1=arrayfun(@minus,int16(I),int16(Ijj));
ch1_=double(ch1);
ch1_1=power(ch1,2);
sum3=gather(sum(ch1_1(:)));
MSE=sum3/pc;
PSNR=10*log10((255^2)/MSE);
index_0=find(ch1(:)~=0);
count=length(index_0);
cwl=count/pc;%错误率
zql=(pc-count)/pc;%正确率
jg(b,1)=count;
jg(b,2)=PSNR;
jg(b,3)=cwl;

% IIJM=gather(Ij);
% IIjiemi=gather(Ijj);
% imwrite(IIJM,'result\goldhill\jiami.png');
% imwrite(IIjiemi,'result\goldhill\jiemi.png');