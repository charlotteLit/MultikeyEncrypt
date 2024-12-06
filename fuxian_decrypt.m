function [recover_I, ext_D] = fuxian_decrypt(encrypt_I, bs, b, row, col)

    rand_num = 2^b;
    rng(b);
    rand_seq = floor(rand(rand_num, row, col) * 255);

    % �ֿ���
    m = floor(row/bs);
    n = floor(col/bs);

    decrypt_I = encrypt_I;
    % �Ȼָ���һ�е�һ��
    key = squeeze(rand_seq(1,:,:));
    decrypt_I(1:bs,:) = bitxor(decrypt_I(1:bs,:), key(1:bs,:));
    decrypt_I(bs+1:row,1:bs) = bitxor(decrypt_I(bs+1:row,1:bs),key(bs+1:row,1:bs));

    % ��������Կ���������ͼ����򣬵õ�rand_num*row*(col)��С�ľ���
    multi_encrypt=permute(repmat(encrypt_I,1,1,rand_num),[3,1,2]);
    xor_array=arrayfun(@bitxor,multi_encrypt,rand_seq);

    % ����ԭ�ĵ�һ�е�һ�и�ֵ�����������
    for i=1:rand_num
        % ��һ�е�һ���м���
        xor_array(i,1:bs,:) = decrypt_I(1:bs,:);
        xor_array(i,bs+1:row,1:bs) = decrypt_I(bs+1:row,1:bs);
    end

    % ����xor_array�е�һ�е�һ���м����ѽ��ܣ�����δ���ܣ�Ѱ����Сf
    sum_block_l = zeros(rand_num,m,n);
    for index = 1:rand_num
        
        for i = 2:m
            b_x = (i-1) * bs+1;
            e_x = bs*i;
            for j = 2:n
                b_y = (j-1) * bs+1;
                e_y = bs*j;
                % ���ڼ���
                % ˮƽ����
                hor_l = zeros();
                cnt1_l = 0;
                for r=b_x:e_x
                    for c=b_y:e_y
                        hor_l(cnt1_l+1)=xor_array(index,r,c)-xor_array(index,r,c-1);
                        cnt1_l = cnt1_l + 1;
                    end
                end
                % ��ֱ����
                ver_l = zeros();
                cnt2_l = 0;
                for r=b_x:e_x
                    for c=b_y:e_y
                        ver_l(cnt2_l+1)=xor_array(index,r,c)-xor_array(index,r-1,c);
                        cnt2_l = cnt2_l + 1;
                    end
                end
                % �������ֵ֮��
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

    % ������Сֵ���õ����ڽ��ܵ�decrypt_mask
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
    % ��ԭ�м�ͼ��
    recover_I = arrayfun(@bitxor,recover_I,decry);
    % ��ԭ��Ե����
    % �Ȼָ���һ�е�һ��
    recover_I(1:bs,:) = decrypt_I(1:bs,:);
    recover_I(bs+1:row,1:bs) = decrypt_I(bs+1:row,1:bs);