function [recover_I, ext_D] = fuxian_decrypt(encrypt_I, bs, b, row, col)

    rand_num = 2^b;
    rng(b);
    rand_seq = floor(rand(rand_num, row, col) * 255);

    % 分块数
    m = floor(row/bs);
    n = floor(col/bs);

    decrypt_I = encrypt_I;
    % 先恢复第一行第一列
    key = squeeze(rand_seq(1,:,:));
    decrypt_I(1:bs,:) = bitxor(decrypt_I(1:bs,:), key(1:bs,:));
    decrypt_I(bs+1:row,1:bs) = bitxor(decrypt_I(bs+1:row,1:bs),key(bs+1:row,1:bs));

    % 将所有密钥依次与加密图像异或，得到rand_num*row*(col)大小的矩阵
    multi_encrypt=permute(repmat(encrypt_I,1,1,rand_num),[3,1,2]);
    xor_array=arrayfun(@bitxor,multi_encrypt,rand_seq);

    % 将复原的第一行第一列赋值到异或数组中
    for i=1:rand_num
        % 第一行第一列中间列
        xor_array(i,1:bs,:) = decrypt_I(1:bs,:);
        xor_array(i,bs+1:row,1:bs) = decrypt_I(bs+1:row,1:bs);
    end

    % 现在xor_array中第一行第一列中间列已解密，其余未解密，寻找最小f
    sum_block_l = zeros(rand_num,m,n);
    for index = 1:rand_num
        
        for i = 2:m
            b_x = (i-1) * bs+1;
            e_x = bs*i;
            for j = 2:n
                b_y = (j-1) * bs+1;
                e_y = bs*j;
                % 块内计算
                % 水平方向
                hor_l = zeros();
                cnt1_l = 0;
                for r=b_x:e_x
                    for c=b_y:e_y
                        hor_l(cnt1_l+1)=xor_array(index,r,c)-xor_array(index,r,c-1);
                        cnt1_l = cnt1_l + 1;
                    end
                end
                % 竖直方向
                ver_l = zeros();
                cnt2_l = 0;
                for r=b_x:e_x
                    for c=b_y:e_y
                        ver_l(cnt2_l+1)=xor_array(index,r,c)-xor_array(index,r-1,c);
                        cnt2_l = cnt2_l + 1;
                    end
                end
                % 计算绝对值之和
                sum_l = 0;
                for cnt = 1:cnt1_l
                    sum_l = sum_l + abs(hor_l(cnt)) + abs(ver_l(cnt));
                end
                sum_block_l(index,i,j) = sum_l;
            end
        end

    end % end of index

    sum_block = sum_block_l;
    [~,index]=min(sum_block,[],1);
    index = squeeze(index);

    % 根据最小值，得到用于解密的decrypt_mask
    decry = zeros(row, col);
    ext_D = zeros();
    ed = 0;
    for i=2:m
        begin_x = (i-1) * bs;
        x_range = begin_x+1:begin_x+bs;
        for j = 2:n
            begin_y = (j-1) * bs;
            y_range = begin_y+1:begin_y+bs;
            ind = index(i,j);
            extract = decimal_binary(ind-1, b);
            ext_D(ed+1:ed+b) = extract(1:b);
            ed = ed + b;
            key = squeeze(rand_seq(ind,:,:));
            decry(x_range,y_range) = key(x_range,y_range);
        end
    end

    recover_I = encrypt_I;
    % 复原中间图像
    recover_I = arrayfun(@bitxor,recover_I,decry);
    % 复原边缘区域
    % 先恢复第一行第一列
    recover_I(1:bs,:) = decrypt_I(1:bs,:);
    recover_I(bs+1:row,1:bs) = decrypt_I(bs+1:row,1:bs);