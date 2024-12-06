function [encrypt_I, emd_data] = encrypt(double_I, bs, b, D, row, col)

    rand_num = 2^b;
    rng(5);
    rand_seq = floor(rand(rand_num, row, col*2) * 255);

    % ������Ϣ���
    d = 0;
    encrypt_I = double_I;

    % �ֿ���
    m = floor(row/bs);
    n = floor(col/bs);

    % ��һ�е�һ�лҼ�
    % ���ɼ�������mask
    mask = zeros(row,col*2);

    % ���ܻҼ���ѡ����Կ1
    key = squeeze(rand_seq(1,:,:));

    % ��һ�е�һ��
    mask(1:bs,:) = key(1:bs,:);
    mask(bs+1:row,1:bs) = key(bs+1:row,1:bs);
    % �м���
    mask(bs+1:row,col+1:col+bs) = key(bs+1:row,col+1:col+bs);


    emd_data = zeros();
    emd = 0;

    % ����mask����
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

    % ���ܰ׼�
    encrypt_I = arrayfun(@bitxor,encrypt_I,mask);