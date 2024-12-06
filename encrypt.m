function [encrypt_I, emd_data] = encrypt(double_I, bs, b, D, row, col)

    rand_num = 2^b;
    rng(5);
    rand_seq = floor(rand(rand_num, row, col*2) * 255);

    % 隐藏信息标记
    d = 0;
    encrypt_I = double_I;

    % 分块数
    m = floor(row/bs);
    n = floor(col/bs);

    % 第一行第一列灰集
    % 生成加密掩码mask
    mask = zeros(row,col*2);

    % 加密灰集，选定密钥1
    key = squeeze(rand_seq(1,:,:));

    % 第一行第一列
    mask(1:bs,:) = key(1:bs,:);
    mask(bs+1:row,1:bs) = key(bs+1:row,1:bs);
    % 中间列
    mask(bs+1:row,col+1:col+bs) = key(bs+1:row,col+1:col+bs);


    emd_data = zeros();
    emd = 0;

    % 完善mask内容
    for i = 2:m
        begin_x = (i-1) * bs;
        x_range = begin_x+1:begin_x+bs;
        for j = 2:n
            ed = binary_decimal(D(d+1:d+b)) + 1;
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