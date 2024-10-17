%% 准备原始图像，分两幅

clc
clear

I = imread('./测试图像/Lena.tiff');
origin_I = double(I);

inter_I = origin_I(1:510,1:510);
[row,col] = size(inter_I);

double_I = [inter_I, inter_I];
figure(1);
imshow(double_I,[]);
title('拼接图像');
imwrite(uint8(double_I),'catch.png');

%% 图像加密

% 分块大小
bs = 3;

% 分块数
m = floor(row/bs);
n = floor(col/bs);

% 嵌入数据生成
num_D = 1000000;
rng(0);
D = round(rand(1,num_D)*1);

% 生成2^b个随机序列，大小为两幅图片水平排列
% 每块嵌入5 bit数据
b = 5;
rand_num = 2^b;
rng(b);
rand_seq = floor(rand(rand_num, row, col*2) * 255);

% 图像加密 & 数据嵌入
% 隐藏信息标记
d = 0;
encrypt_I = double_I;

% 第一行第一列灰集
% 生成加密掩码mask
mask = zeros(row,col*2);

% 加密灰集，选定密钥1
key = squeeze(rand_seq(1,:,:));

% 第一行第一列
mask(1:bs,:) = key(1:bs,:);
mask(bs+1:row,1:bs) = key(bs+1:row,1:bs);
% 中间列
mask(bs+1:row,col+1:col+bs) = key(bs+1:510,col+1:col+bs);

% 最后列，最后行
% mask(bs+1:510,end-2:end) = key(bs+1:510,end-2:end);
% mask(510:512,bs+1:end) = key(510:512,bs+1:end);

emd_data = zeros();
emd = 0;

% 完善mask内容
for i = 2:m
    begin_x = (i-1) * bs;
    x_range = begin_x+1:begin_x+bs;
    for j = 2:n
        ed = binary_decimal([0, 0, 0, D(d+1:d+b)]) + 1;
        emd_data(emd+1:emd+b) = D(d+1:d+b);
        d = d + b;
        emd = emd + b;
        key = squeeze(rand_seq(ed, :, :));

        begin_y = (j-1) * bs;
        begin_y1 = begin_y + col;
        y_range = begin_y+1:begin_y+bs;
        y1_range = begin_y1+1:begin_y1+bs;
        mask(x_range,y_range) = key(x_range,y_range);
        mask(x_range, y1_range) = key(x_range, y1_range);
    end  
end

% 加密白集
encrypt_I = arrayfun(@bitxor,encrypt_I,mask);

figure(2);
imshow(encrypt_I,[]);
title('加密图像2');
imwrite(uint8(encrypt_I),'encrypt.png');

%% 数据提取 & 复原图像

decrypt_I = encrypt_I;
% 先恢复第一行第一列
key = squeeze(rand_seq(1,:,:));
decrypt_I(1:bs,:) = bitxor(decrypt_I(1:bs,:), key(1:bs,:));
decrypt_I(bs+1:row,1:bs) = bitxor(decrypt_I(bs+1:row,1:bs),key(bs+1:row,1:bs));
% 中间列
decrypt_I(bs+1:row,col+1:col+bs) = bitxor(decrypt_I(bs+1:row,col+1:col+bs),key(bs+1:row,col+1:col+bs));
% 最后列，最后行
% decrypt_I(bs+1:510,end-2:end) = bitxor(decrypt_I(bs+1:510,end-2:end),key(bs+1:510,end-2:end));
% decrypt_I(510:512,bs+1:end) = bitxor(decrypt_I(510:512,bs+1:end),key(510:512,bs+1:end));

% 将所有密钥依次与加密图像异或，得到rand_num*row*(col*2)大小的矩阵
multi_encrypt=permute(repmat(encrypt_I,1,1,rand_num),[3,1,2]);
xor_array=arrayfun(@bitxor,multi_encrypt,rand_seq);

% 将复原的第一行第一列赋值到异或数组中
for i=1:rand_num
    % 第一行第一列中间列
    xor_array(i,1:bs,:) = decrypt_I(1:bs,:);
    xor_array(i,bs+1:row,1:bs) = decrypt_I(bs+1:row,1:bs);
    xor_array(i,bs+1:row,col+1:col+bs) = decrypt_I(bs+1:row,col+1:col+bs);
end

% 现在xor_array中第一行第一列中间列已解密，其余未解密，寻找最小f
sum_block_l = zeros(rand_num,m,n);
sum_block_r = zeros(rand_num,m,n);
for index = 1:rand_num
    
    for i = 2:m
        b_x = (i-1) * bs+1;
        e_x = bs*i;
        for j = 2:n
            b_y = (j-1) * bs+1;
            e_y = bs*j;
            b_y1 = b_y + col;
            e_y1 = e_y + col;
            % 块内计算
            % 水平方向左图
            hor_l = zeros();
            cnt1_l = 0;
            for r=b_x:e_x
                for c=b_y:e_y
                    hor_l(cnt1_l+1)=xor_array(index,r,c)-xor_array(index,r,c-1);
                    cnt1_l = cnt1_l + 1;
                end
            end
            % 竖直方向左图
            ver_l = zeros();
            cnt2_l = 0;
            for r=b_x:e_x
                for c=b_y:e_y
                    ver_l(cnt2_l+1)=xor_array(index,r,c)-xor_array(index,r-1,c);
                    cnt2_l = cnt2_l + 1;
                end
            end
            % 水平方向右图
            hor_r = zeros();
            cnt1_r = 0;
            for r=b_x:e_x
                for c=b_y1:e_y1
                    hor_r(cnt1_r+1)=xor_array(index,r,c)-xor_array(index,r,c-1);
                    cnt1_r = cnt1_r + 1;
                end
            end
            % 竖直方向右图
            ver_r = zeros();
            cnt2_r = 0;
            for r=b_x:e_x
                for c=b_y1:e_y1
                    ver_r(cnt2_r+1)=xor_array(index,r,c)-xor_array(index,r-1,c);
                    cnt2_r = cnt2_r + 1;
                end
            end
            % 计算绝对值之和
            % 左图
            sum_l = 0;
            for cnt = 1:cnt1_l
                sum_l = sum_l + abs(hor_l(cnt)) + abs(ver_l(cnt));
            end
            sum_block_l(index,i,j) = sum_l;
            % 右图
            sum_r = 0;
            for cnt = 1:cnt1_r
                sum_r = sum_r + abs(hor_r(cnt)) + abs(ver_r(cnt));
            end
            sum_block_r(index,i,j) = sum_r;
        end
    end
end % end of index

% 两组对照相加取最小值对应块
sum_block = sum_block_l + sum_block_r;
[value,index]=min(sum_block,[],1);
index = squeeze(index);

% 根据最小值，得到用于解密的decrypt_mask
decry = zeros(row, col*2);
for i=2:m
    begin_x = (i-1) * bs;
    x_range = begin_x+1:begin_x+bs;
    for j = 2:n
        begin_y = (j-1) * bs;
        begin_y1 = begin_y + col;
        y_range = begin_y+1:begin_y+bs;
        y1_range = begin_y1+1:begin_y1+bs; 
        ind = index(i,j);
        key = squeeze(rand_seq(ind,:,:));
        decry(x_range,y_range) = key(x_range,y_range);
        decry(x_range,y1_range)= key(x_range,y1_range);
    end
end


recover_I = encrypt_I;
% 复原中间图像
recover_I = arrayfun(@bitxor,recover_I,decry);
% 复原边缘区域
% 先恢复第一行第一列
key = squeeze(rand_seq(1,:,:));
recover_I(1:bs,:) = decrypt_I(1:bs,:);
recover_I(bs+1:row,1:bs) = decrypt_I(bs+1:row,1:bs);
% 中间列
recover_I(bs+1:row,col+1:col+bs) = decrypt_I(bs+1:row,col+1:col+bs);
% 最后列，最后行
% recover_I(bs+1:510,end-2:end) = decrypt_I(bs+1:510,end-2:end);
% recover_I(510:512,bs+1:end) = decrypt_I(510:512,bs+1:end);

figure(3);
imshow(recover_I,[]);
title('解密图像');
recover_single = recover_I(:,1:col);
imwrite(uint8(recover_single),'recover.png');

%% 计算MSE均方误差，错误率，正确率
pixel = row*col;
dif=arrayfun(@minus,int16(inter_I),int16(recover_single));
d_dif=double(dif);
squ_dif=power(d_dif,2);
sum_dif=sum(squ_dif(:));
MSE=sum_dif/pixel;
PSNR=10*log10((255^2)/MSE);

index_error=find(dif(:)~=0);
count=length(index_error);
error_rate=count/pixel;% 错误率
correct_rate=(pixel-count)/pixel;% 正确率
