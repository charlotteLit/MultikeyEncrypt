function [MSE, PSNR, ER, BPP] = caulate(row, col, origin_I, recover_I, emd_D)
    % ����MSE�����������ʣ���ȷ��
    pixel = row*col;
    dif=arrayfun(@minus,int16(origin_I),int16(recover_I));
    d_dif=double(dif);
    squ_dif=power(d_dif,2);
    sum_dif=sum(squ_dif(:));

    MSE=sum_dif/pixel;
    PSNR=10*log10((255^2)/MSE);

    index_error=find(dif(:)~=0);
    count=length(index_error);
    ER=count/pixel; % ������
    BPP = emd_D/pixel; % Ƕ����