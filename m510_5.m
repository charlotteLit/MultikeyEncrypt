I=imread('Barbara.tif');
%I=rgb2gray(I);
I=I(1:510,1:510);
[r1,c1]=size(I);
m=fix(r1/5);
n=fix(c1/5);
I=gpuArray(I);
 b=5;
%imshow(I);
%����ҪǶ�������ͼ�񣬺ڰ�ͼ��
J1=imread('hei1.jpg');
J1=imresize(J1,[b,10201]);
J2=rgb2gray(J1);
J=im2bw(J2);
[r2,c2]=size(J);
J=im2uint8(J);%����ff����������logical����
J=imnoise(J,'gaussian',0.02);
J=imnoise(J,'poisson');
J=im2bw(J);
% J=gpuArray(J);


jg=zeros(11,3);
...........����wd���������............

%b=11;
wd=2^b;
rng('shuffle');
p=floor(255*rand([wd,510,510]));


..........��ȷ��ÿ��ҪǶ��ı�������Ȼ��ѭ������Ƕ�룬�ڴ�ѭ����߼���ÿ�ε���Կֵ��Ȼ�������Կֵ�Ĳ�ͬ��ѡ��ͬ��������������......
zc=m*n*b;    %����ͼ������ҪǶ��ı�����
mc=r2*c2;    %����ͼ��ӵ�еı�����



tic;

    %�ȰѴ�Ƕ�������ͼ�����δ���һ��һά����
    JJ=zeros(1,m*n*b);
    %�˴�ת��
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
    %����ÿһ���Ӧ��y_ֵ����505*505��С���������
    mask=zeros(510,510);
    maxindex=length(y_);
    
    GetMask5(maxindex,p,y_,mask,wd);
%     for index=1:maxindex;
%         ȡ����һ��������Ϣ��Ӧ��ʮ������
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


% �ÿ���һ�к͵�һ�н��м���
  mmask=mask(6:510,6:510);
  Ip=I(6:510,6:510);

Ijp=arrayfun(@bitxor,Ip,mmask);%�������6:510,6:510����

Ij=I;
Ij(6:510,6:510)=Ijp;
% 

% �Ե�һ�е�һ��ʹ��ѡ������Կ���м���
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


subplot(221);imshow(I);title('ԭʼͼ��');
subplot(222);imshow(Ij);title('����ͼ��');


toc;
...............�Ѽ��ܺ��ͼ��ֱ���ÿһ�������������򣬲����ÿ�ε����ھ���֮��...................    

tic;

%������֪key���ܵ�һ�е�һ��
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




%1024��������������ͼ��������

Ij1=permute(repmat(Ij,1,1,wd),[3,1,2]);  %�Ѽ��ܺ�ͼ��Ij����1024�� �����г�1024*510*510
Ij_=arrayfun(@bitxor,Ij1,uint8(p));      %����ͼ������Կ�������1024*510*510ͼƬ
Ij__=gather(Ij_);
%Ij__=Ij__(:,6:510,6:510);
Ijjrc=permute(repmat(Ijj,1,1,wd),[3,1,2]);%�������һ�к͵�һ�к��Ijj����1024��

Ijj_=gather(Ijj);
% for i=1:1024;
%     i
%     Ijrep(i,6:510,6:510)=Ij__(i,:,:);
%     Ijrep(i,1:5,1:5)=Ijj_(1:5,1:5);
% end

%IjjΪ���ֽ��ܵļ���ͼ��
%imshow(Ijj_);
Ijrep=permute(repmat(Ijj,1,1,wd),[3,1,2]);%IjrepΪ��һ�е�һ�н������ʣ�ಿ��Ϊ����ͼ���1024��ͼ��
% %��Ijrep��1024��6:510,6:510�滻��Ijbit_��6:510,6:510
for i=1:2^b
    i
    Ijrep(i,6:510,6:510)=Ij__(i,6:510,6:510);
end
% 
% Ijbit_=gather(Ijbit);
% Ijj_=gather(Ijj);
% 
Ijrep_=gather(Ijrep);
%��ʱIjrep_����1024��ͼ ÿ��ͼ��һ�е�һ��Ϊ������ɲ��� 6:510,6:510Ϊ����ͼ������Կ������ͼ��

% Ijrep(1:5,:)=Ijj_(1:5,:);
% imshow(Ijrep);





%------Ijrep_Ϊ�������(1:5,1:5)�����󲿷�(6:510,6:510)ƴ�ɵ�1024*510*510ͼ��

% Ijjrep=repmat(Ijj,1,1,1024);
% Ijjrep=permute(Ijjrep,[3,1,2]);%Ijjrep�д洢1024����õ�һ�е�һ�в���ʣ�ಿ��Ϊ����ͼ�ε�ͼƬ
% Ijjrep_=gather(Ijjrep);


%��ʱIj__�д洢�ż��ܺ�ͼ��Ij�ģ�6:510,6:510����wd��p�ģ�6:510,6:510������Ľ��

%���������Ӧ�����Ѿ������Ѿ���õĵ�һ�е�һ�еĿ�

sum1=zeros(wd,102,102);
for index_=1:wd
    
   sum0=zeros(102,102);
   for i=2:102
       
       for j=2:102
           
            x1=5*(i-1)+1;
            x5=5*i;
            y1=5*(j-1)+1;
            y5=5*j;
            %һ��5*5�Ŀ�ֻ�ܼ���4*5=20��
            d1=zeros(30,1);
            count1=1;
            %�õ�һ�еĿ����ο� ���׼ȷ��
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
%sum1�д洢��1024��102*102������
%ȥ��sum1�е�һ�е�һ��

toc
%     fun=@(block_struct) abs(block_struct.data(1,2)-block_struct.data(1,1))+abs(block_struct.data(1,3)-block_struct.data(1,2))+abs(block_struct.data(2,2)-block_struct.data(2,1))+abs(block_struct.data(2,3)-block_struct.data(2,2))+abs(block_struct.data(3,2)-block_struct.data(3,1))+abs(block_struct.data(3,3)-block_struct.data(3,2))+abs(block_struct.data(2,1)-block_struct.data(1,1))+abs(block_struct.data(2,2)-block_struct.data(1,2))+abs(block_struct.data(2,3)-block_struct.data(1,3))+abs(block_struct.data(3,1)-block_struct.data(2,1))+abs(block_struct.data(3,2)-block_struct.data(2,2))+abs(block_struct.data(3,3)-block_struct.data(2,3));
%     result=blockproc(squeeze(Ij_(i,:,:)),[3,3],fun);
%     sum1(i,:,:)=result;


%�ֱ��ҳ�ÿ�����Ӧ����С�ͣ�������������ͼ��Ķ�Ӧ�������򣬴Ӷ������ȷ��ԭͼ��

%sum1�д洢1024��102*102 ������


 %[Y,U]=max(sum1)%����������Y��U��Y������¼A��ÿ�е����ֵ��U������¼ÿ�����ֵ���кš�
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

%��101*101����Կ�������һ��505*505��mask_


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
%Ijq=Ij(6:510,6:510);%ȡ����ͼ���ʣ�ಿ�֣���ȥ��һ�� ��һ�У�
%�ҳ�sum11��ÿһ���Ӧ����Կ���
[qq,mk]=min(sum1,[],1);
mk=reshape(mk,102,102);
%�˴�ת��
mk=mk';
mask_=zeros(510,510);
GetMask5(maxindex,p,mk,mask_,wd);
mmask_=mask_(6:510,6:510);
% mk=mk';
% mk=mk(:);%��mk�ų�һά��������
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

%����ͼ����ƴ�ɵ�mask_���
Ijjp=arrayfun(@bitxor,Ijp,mmask_);
Ijj(6:510,6:510)=Ijjp;
subplot(223);imshow(Ijj);title('����ͼ��');
toc

%������ܺ��ͼ��Ĵ�����
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
cwl=count/pc;%������
zql=(pc-count)/pc;%��ȷ��
jg(b,1)=count;
jg(b,2)=PSNR;
jg(b,3)=cwl;

% IIJM=gather(Ij);
% IIjiemi=gather(Ijj);
% imwrite(IIJM,'result\goldhill\jiami.png');
% imwrite(IIjiemi,'result\goldhill\jiemi.png');