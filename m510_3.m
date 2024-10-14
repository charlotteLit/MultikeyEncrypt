I=imread('Goldhill.png');
% ��ɫͼ��ת�Ҷ�ͼ��
% I=rgb2gray(I);
I=I(1:510,1:510);
[r1,c1]=size(I);
m=fix(r1/3);
n=fix(c1/3);
% ��Ӳ�������
% I=imnoise(I,'poisson');

% תΪGPU���飬�������� GPU ���ٴ���
I=gpuArray(I);
b=5;

%����ҪǶ�������ͼ�񣬺ڰ�ͼ��
J1=imread('hei1.jpg');
J2=rgb2gray(J1);

% ������ֵתΪ��ֵͼ��
J=imbinarize(J2);
[r2,c2]=size(J);
J=im2uint8(J);%����ff����������logical����
% ��Ӹ�˹������0.02������ǿ��
J=imnoise(J,'gaussian',0.02);
% ��������
J=imnoise(J,'poisson');
J=imbinarize(J);
J=J';
% J=gpuArray(J);


jg=zeros(11,3);
...........����wd���������............

%b=11;
wd=2^b;
% �����������
rng(b,'v4');
p=floor(255*rand([wd,510,510]));


...��ȷ��ÿ��ҪǶ��ı�������Ȼ��ѭ������Ƕ�룬�ڴ�ѭ����߼���ÿ�ε���Կֵ
...Ȼ�������Կֵ�Ĳ�ͬ��ѡ��ͬ��������������...
zc=m*n*b;    %����ͼ������ҪǶ��ı�����
mc=r2*c2;    %����ͼ��ӵ�еı�����


% ����������
tic;

    %�ȰѴ�Ƕ�������ͼ�����δ���һ��һά����
    JJ=zeros(1,m*n*b);
    %�˴�ת��
    J=J';
    JJ(1,1:r2*c2)=J(:);
    JJ=reshape(JJ,b,[]);
    [j_c,j_r]=size(JJ);
    x=zeros(b,1);
    for i=1:b 
        x(i,1)=2^(i-1);
    end
    % ����x����һ��b*j_r��С������
    x=repmat(x,1,j_r);
    % �����Ȩ��
    y=arrayfun(@mtimes,JJ,x);
    y=num2cell(y,1);
    % �����뵥Ԫ��������Ԫ�����
    F=@(x)sum(x{:});
    % ����y��ÿ����Ԫ�ĺ�
    y_=arrayfun(F,y);
    %y_=y_+1;
    %����ÿһ���Ӧ��y_ֵ����505*505��С���������
    mask=zeros(510,510);
    maxindex=length(y_);
    yy=y_+1;
    yy=reshape(yy,170,170);
    yy=yy';
    for index=1:maxindex
        num_=y_(1,index)+1;
        q=squeeze(p(num_,:,:));
        r_index=ceil(index/170)*3-2;
        c_index_=mod(index-170*fix(index/170),170);
        if c_index_==0
            c_index_=170;
        end
        c_index=c_index_*3-2;
        mask(r_index:r_index+2,c_index:c_index+2)=q(r_index:r_index+2,c_index:c_index+2);
    end
    
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
  mmask=mask(4:510,4:510);
  Ip=I(4:510,4:510);

Ijp=arrayfun(@bitxor,Ip,mmask);%�������4:510,4:510����

Ij=I;
Ij(4:510,4:510)=Ijp;
% 

% �Ե�һ�е�һ��ʹ��ѡ������Կ���м���
key1=p(1,:,:); 
key1=uint8(key1);
key1=reshape(key1,510,510);
for i=1:170
   i1=(i-1)*3+1;
   i2=i*3;
   Ij(1:3,i1:i2)=bitxor(key1(1:3,i1:i2),I(1:3,i1:i2));
end

key2=p(2,:,:); 
key2=uint8(key2);
key2=reshape(key2,510,510);
for i=2:170
   i1=(i-1)*3+1;
   i2=i*3;
   Ij(i1:i2,1:3)=bitxor(key2(i1:i2,1:3),I(i1:i2,1:3));
end


subplot(221);imshow(I);title('ԭʼͼ��');
subplot(222);imshow(Ij);title('����ͼ��');


toc;
...............�Ѽ��ܺ��ͼ��ֱ���ÿһ�������������򣬲����ÿ�ε����ھ���֮��...................    

tic;

%������֪key���ܵ�һ�е�һ��
Ijj=Ij;

for i=1:170
   i1=(i-1)*3+1;
   i2=i*3;
   Ij1r(1:3,i1:i2)=bitxor(key1(1:3,i1:i2),Ij(1:3,i1:i2));
end
Ijj(1:3,:)=Ij1r(1:3,:);


for i=2:170
    i1=(i-1)*3+1;
    i2=i*3;
    Ij1r(i1:i2,1:3)=bitxor(key2(i1:i2,1:3),Ij(i1:i2,1:3));
end
Ijj(:,1:3)=Ij1r(:,1:3);


%subplot(223);imshow(Ijj);title('����ͼ��');


%1024��������������ͼ��������

Ij1=permute(repmat(Ij,1,1,wd),[3,1,2]);  %�Ѽ��ܺ�ͼ��Ij����1024�� �����г�1024*510*510
Ij_=arrayfun(@bitxor,Ij1,uint8(p));      %����ͼ������Կ�������1024*510*510ͼƬ
Ij__=gather(Ij_);
%Ij__=Ij__(6:510,6:510);
Ijjrc=permute(repmat(Ijj,1,1,wd),[3,1,2]);%�������һ�к͵�һ�к��Ijj����1024��

Ijj_=gather(Ijj);
% for i=1:1024;
%     Ijrep(i,6:510,6:510)=Ij__(i,6:510,6:510);
%     Ijrep(i,1:5,1:5)=Ijj_(1:5,1:5);
% end
%IjjΪ���ֽ��ܵļ���ͼ��
%imshow(Ijj_);
Ijrep=permute(repmat(Ijj,1,1,wd),[3,1,2]);%IjrepΪ��һ�е�һ�н������ʣ�ಿ��Ϊ����ͼ���1024��ͼ��
% %��Ijrep��1024��6:510,6:510�滻��Ijbit_��6:510,6:510
for i=1:2^b
    Ijrep(i,4:510,4:510)=Ij__(i,4:510,4:510);
end
Ijrep_=gather(Ijrep);
% 
% Ijbit_=gather(Ijbit);
% Ijj_=gather(Ijj);
% 
% Ijrep_=gather(Ijrep);
%��ʱIjrep����1024��ͼ ÿ��ͼ��һ�е�һ��Ϊ������ɲ��� 6:510,6:510Ϊ����ͼ������Կ������ͼ��

% Ijrep(1:5,:)=Ijj_(1:5,:);
% imshow(Ijrep);





%------Ijrep_Ϊ�������(1:5,1:5)�����󲿷�(6:510,6:510)ƴ�ɵ�1024*510*510ͼ��

% Ijjrep=repmat(Ijj,1,1,1024);
% Ijjrep=permute(Ijjrep,[3,1,2]);%Ijjrep�д洢1024����õ�һ�е�һ�в���ʣ�ಿ��Ϊ����ͼ�ε�ͼƬ
% Ijjrep_=gather(Ijjrep);


%���������Ӧ�����Ѿ������Ѿ���õĵ�һ�е�һ�еĿ�

sum1=zeros(wd,170,170);
sum2=zeros(wd,170,170);
for index_=1:wd
    index_
   sum0=zeros(170,170);
   sum00=zeros(170,170);
   %������
   for i=2:170
       for j=2:170
           
            x1=3*(i-1)+1;
            x3=3*i;
            y1=3*(j-1)+1;
            y3=3*j;
            
            d1=zeros(20,1);
            count1=1;
            
            for i1=x1:x3
                for j1=y1:(y3-1)
                    d1(count1)=double(int16(Ijrep_(index_,i1,j1+1))-int16(Ijrep_(index_,i1,j1)));
                    count1=count1+1;
                end
            end
            
            d2=zeros(20,1);
            count2=1;
            for j2=y1-1:y3
                for i2=x1-1:(x3-1)
                    d2(count2)=double(int16(Ijrep_(index_,i2+1,j2))-int16(Ijrep_(index_,i2,j2)));
                    count2=count2+1;
                end
            end
            
            s=0;
            for i3=1:20
                s=s+abs(d1(i3))+abs(d2(i3));
            end
            sum0(i,j)=s;
        end
    end
   sum1(index_,:,:)=sum0;
  
   %������
   for i=170:-1:2
       for j=170:-1:2
           
            x1=3*(i-1)+1;
            x3=3*i;
            y1=3*(j-1)+1;
            y3=3*j;
            
            d3=zeros(20,1);
            count1=1;
            
            for i1=x1:x3
                for j1=y1:(y3-1)
                    d3(count1)=double(int16(Ijrep_(index_,i1,j1+1))-int16(Ijrep_(index_,i1,j1)));
                    count1=count1+1;
                end
            end
            
            d4=zeros(20,1);
            count2=1;
            for j2=y1-1:y3
                for i2=x1-1:(x3-1)
                    d4(count2)=double(int16(Ijrep_(index_,i2+1,j2))-int16(Ijrep_(index_,i2,j2)));
                    count2=count2+1;
                end
            end
            
            ss=0;
            for i3=1:20
                ss=ss+abs(d3(i3))+abs(d4(i3));
            end
            sum00(i,j)=ss;
        end
    end
   sum2(index_,:,:)=sum0;
end
%sum1�д洢��1024��102*102������
%ȥ��sum1�е�һ�е�һ��
sum=sum1+sum2;
toc
%     fun=@(block_struct) 
%   abs(block_struct.data(1,2)-block_struct.data(1,1))+abs(block_struct.data(1,3)-block_struct.data(1,2))+abs(block_struct.data(2,2)-block_struct.data(2,1))+abs(block_struct.data(2,3)-block_struct.data(2,2))+abs(block_struct.data(3,2)-block_struct.data(3,1))+abs(block_struct.data(3,3)-block_struct.data(3,2))+abs(block_struct.data(2,1)-block_struct.data(1,1))+abs(block_struct.data(2,2)-block_struct.data(1,2))+abs(block_struct.data(2,3)-block_struct.data(1,3))+abs(block_struct.data(3,1)-block_struct.data(2,1))+abs(block_struct.data(3,2)-block_struct.data(2,2))+abs(block_struct.data(3,3)-block_struct.data(2,3));
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
[qq,mk]=min(sum,[],1);
mk=reshape(mk,170,170);
%  mk(3,14)=513;
%  mk(5,5)=912;
%�˴�ת��
%mk=mk';
mask_=zeros(510,510);
 for index__=1:n*m
        r_index_=ceil(index__/170);
        c_index_=mod(index__-170*fix(index__/170),170);
        if c_index_==0
            c_index_=170;
        end
        r_index=r_index_*3-2;
        c_index=c_index_*3-2;
        num_=mk(r_index_,c_index_);
        q=squeeze(p(num_,:,:));
        mask_(r_index:r_index+2,c_index:c_index+2)=q(r_index:r_index+2,c_index:c_index+2);
 end
mmask_=mask_(4:510,4:510);

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
Ijj(4:510,4:510)=Ijjp;
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
%PSNR=10*log10((255^2)/MSE);
index_0=find(ch1(:)~=0);
count=length(index_0);
cwl=count/pc;%������
zql=(pc-count)/pc;%��ȷ��
% jg(b,1)=count;
% jg(b,2)=PSNR;
% jg(b,3)=cwl;

IIjiemi=gather(Ijj);
imwrite(IIjiemi,'0226goldhill_3_5.png');

% yy_=yy(2:170,2:170);
% mk_=mk(2:170,2:170);
% find(yy_~=mk_);
% [mrows,mcols]=find(yy_~=mk_);

% cc=0;
% for i=1:56
%     if mrows(i)<11
%         cc=cc+1;
%     end
% 
% end
% cc
% 
% for i=1:56
%     if mrows(i)==2;
%         break
%     else i=6219;
%     end
%     
% end
% i
% 
% IIJM=gather(Ij);
% IIjiemi=gather(Ijj);
% %imwrite(IIJM,'jiami10016.png');
% imwrite(IIjiemi,'ct1jiemi.png');
