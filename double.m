I = imread('./����ͼ��/Lena.tiff');
origin_I = double(I);
[row,col] = size(origin_I);

double_I = [origin_I, origin_I];
figure(1);
imshow(double_I,[]);
title('ƴ��ͼ��');

% �ֿ��С
bs = 3;

% �ֿ���
m = floor(row/bs);
n = floor(col/bs);

% Ƕ����������
num_D = 1000000;
rng(0);
D = round(rand(1,num_D)*1);

% ����2^b��������У���СΪ����ͼƬˮƽ����
% ÿ��Ƕ��5 bit����
b = 5;
rand_num = 2^b;
rng(b);
rand_seq = floor(rand(rand_num, row, col*2) * 255);

% ͼ����� & ����Ƕ��
% ������Ϣ���
d = 0;
encrypt_I = double_I;

% ��һ�е�һ�лҼ�
% ���ɼ�������mask
mask = zeros(row,col*2);

% ���ܻҼ���ѡ����Կ1
key = squeeze(rand_seq(1,:,:));

% ��һ�е�һ��
mask(1:bs,:) = key(1:bs,:);
mask(bs+1:end,1:bs) = key(bs+1:end,1:bs);
% �м��У�����У������
mask(bs+1:510,511:515) = key(bs+1:510,511:515);
mask(bs+1:510,end-2:end) = key(bs+1:510,end-2:end);
mask(510:512,bs+1:end) = key(510:512,bs+1:end);

emd_data = zeros();
emd = 0;

% ����mask����
for i = 2:m
    begin_x = (i-1) * bs;
    x_range = begin_x+1:begin_x+bs;
    for j = 2:n
        ed = bin2dec([0, 0, 0, D(d+1:d+5)]) + 1;
        d = d + 5;
        emd_data(emd+1:emd+5) = D(d+1:d+5);
        emd = emd + 5;
        key = squeeze(rand_seq(ed, :, :));

        begin_y = (j-1) * bs;
        begin_y1 = begin_y + col;
        y_range = begin_y+1:begin_y+bs;
        y1_range = begin_y1+1:begin_y1+bs;
        mask(x_range,y_range) = key(x_range,y_range);
        mask(x_range, y1_range) = key(x_range, y1_range);
    end  
end

% ���ܰ׼�
encrypt_I = arrayfun(@bitxor,encrypt_I,mask);

figure(2);
imshow(encrypt_I,[]);
title('����ͼ��2');

% ������ȡ & ��ԭͼ��

decrypt_I = encrypt_I;
% �Ȼָ���һ�е�һ��
key = squeeze(rand_seq(1,:,:));
decrypt_I(1:bs,:) = bitxor(decrypt_I(1:bs,:), key(1:bs,:));
decrypt_I(bs+1:510,1:bs) = bitxor(decrypt_I(bs+1:510,1:bs),key(bs+1:510,1:bs));
% �м��У�����У������
decrypt_I(bs+1:510,511:515) = bitxor(decrypt_I(bs+1:510,511:515),key(bs+1:510,511:515));
decrypt_I(bs+1:510,end-2:end) = bitxor(decrypt_I(bs+1:510,end-2:end),key(bs+1:510,end-2:end));
decrypt_I(510:512,bs+1:end) = bitxor(decrypt_I(510:512,bs+1:end),key(510:512,bs+1:end));

% ��������Կ���������ͼ����򣬵õ�rand_num*row*(col*2)��С�ľ���
multi_encrypt=permute(repmat(encrypt_I,1,1,rand_num),[3,1,2]);
xor_array=arrayfun(@bitxor,multi_encrypt,rand_seq);

% ����ԭ�ĵ�һ�е�һ�и�ֵ�����������
for i=1:rand_num
    % ��һ�е�һ���м���
    xor_array(i,1:bs,:) = decrypt_I(1:bs,:);
    xor_array(i,bs+1:end,1:bs) = decrypt_I(bs+1:end,1:bs);
    xor_array(i,bs+1:510,511:515) = decrypt_I(bs+1:510,511:515);
end

% ����xor_array�е�һ�е�һ���ѽ��ܣ�����δ���ܣ�Ѱ����Сf
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
            % ���ڼ���
            % ˮƽ������ͼ
            hor_l = zeros();
            cnt1_l = 0;
            for r=b_x:e_x
                for c=b_y:e_y
                    hor_l(cnt1_l+1)=xor_array(index,r,c)-xor_array(index,r,c-1);
                    cnt1_l = cnt1_l + 1;
                end
            end
            % ��ֱ������ͼ
            ver_l = zeros();
            cnt2_l = 0;
            for r=b_x:e_x
                for c=b_y:e_y
                    ver_l(cnt2_l+1)=xor_array(index,r,c)-xor_array(index,r-1,c);
                    cnt2_l = cnt2_l + 1;
                end
            end
            % ˮƽ������ͼ
            hor_r = zeros();
            cnt1_r = 0;
            for r=b_x:e_x
                for c=b_y1:e_y1
                    hor_r(cnt1_r+1)=xor_array(index,r,c)-xor_array(index,r,c-1);
                    cnt1_r = cnt1_r + 1;
                end
            end
            % ��ֱ������ͼ
            ver_r = zeros();
            cnt2_r = 0;
            for r=b_x:e_x
                for c=b_y1:e_y1
                    ver_r(cnt2_r+1)=xor_array(index,r,c)-xor_array(index,r-1,c);
                    cnt2_r = cnt2_r + 1;
                end
            end
            % �������ֵ֮��
            % ��ͼ
            sum_l = 0;
            for cnt = 1:cnt1_l
                sum_l = sum_l + abs(hor_l(cnt)) + abs(ver_l(cnt));
            end
            sum_block_l(index,i,j) = sum_l;
            % ��ͼ
            sum_r = 0;
            for cnt = 1:cnt1_r
                sum_r = sum_r + abs(hor_r(cnt)) + abs(ver_r(cnt));
            end
            sum_block_r(index,i,j) = sum_r;
        end
    end
end % end of index

% ����������ȡ��Сֵ��Ӧ��
sum_block = sum_block_l + sum_block_r;
[value,index]=min(sum_block,[],1);




