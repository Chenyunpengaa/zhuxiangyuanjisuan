[a,R]=geotiffread('E:\test\npp_ccy\occipNPP_2003_01_01.tif');%�ȵ���ͶӰ��Ϣ��ĳ��Ӱ���·�����У������������������е�һ����
info=geotiffinfo('E:\test\npp_ccy\occipNPP_2003_01_01.tif');%ͬ��
[m,n]=size(a);
warning('off');
nppsum=zeros(m*n,18);%�˴�Ҫ�޸ģ����������д���٣���������8���

for year=2003:2020
    filename=strcat('E:\test\npp_ccy\occipNPP_',int2str(year),'_01_01.tif');%�˴�Ҫ�޸ģ��������ǰ����ÿ�����ֲ�����Ƕȵ����ݣ�ע������ļ����֡�
    data=importdata(filename);
    data=reshape(data,m*n,1);
    nppsum(:,year-2002)=data;%�˴���Ҫ�޸ģ��ҵ������Ǵ�2013��ʼ���˴���Ϊ2012.
end

wcsum=zeros(m*n,18);
for year=2003:2020
    filename=strcat('E:\test\pre_ccy\ocuqi_pre_',int2str(year),'_1.tif');%�˴�Ҫ�޸ģ��������ǰ����ÿ������ر��¶ȵ����ݣ�ע������ļ����֡�
    data=importdata(filename);
    data(m+1,:)=[];
    data(:,n+1)=[];
    data=reshape(data,m*n,1);
    wcsum(:,year-2002)=data;
end


for year=2003:2020
    filename=strcat('E:\test\wendu_units\ocuqi_temp_',int2str(year),'_1.tif');%�˴�Ҫ�޸ģ��������ǰ����ÿ������ر��¶ȵ����ݣ�ע������ļ����֡�
    data=importdata(filename);
    data(m+1,:)=[];
    data(:,n+1)=[];
    data=reshape(data,m*n,1);
    wcsumz(:,year-2002)=data;
end




count = 2020-2003+1;
%����Ժ�������
npp_wc_xgx=zeros(m,n);
npp_wc_p=zeros(m,n);
for i=1:length(nppsum)
    npp=nppsum(i,:);
    if min(npp)>=0 %ע�������NPP����Ч��Χ�Ǵ���0������Լ���������Ч��Χ��С��0�Ļ�������Բ��ü����
        wc=wcsum(i,:);
        wz=wcsumz(i,:);
         %[r2,p2]=partialcorri(npp,wc,wz);
         A1 = npp;
         B1 = wc;
         C1 = wz;
         [rab2,pab2]=corrcoef(A1,B1);
         [rac2,pac2]=corrcoef(A1,C1);
         [rbc2,pbc2]=corrcoef(B1,C1);
         r_ab_c=(rab2(2)-rac2(2)*rbc2(2))/(((1-(rac2(2)^2))^(1/2))*((1-rbc2(2)^2)^(1/2)));
         t = (r_ab_c*sqrt(count-3))/(sqrt(1-(r_ab_c^2)));
         
         npp_wc_xgx(i)=r_ab_c;
         npp_wc_p(i)=t;
    end
end

filename5='E:\test\jiangyu_units\ocu\shiyan_pre5p.tif';%�˴�Ҫ�޸ģ������·��������
filename6='E:\test\jiangyu_units\ocu\������_pre25.tif';%ͬ��


%%���ͼ��
geotiffwrite(filename5,npp_wc_xgx,R,'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag);
geotiffwrite(filename6,npp_wc_p,R,'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag);



npp_wc_xgx(npp_wc_p<2.10092)=NaN;
name1='E:\test\jiangyu_units\ocu\ͨ��������0.05�������ƫ���03-20����ֵ05.tif';
geotiffwrite(name1,npp_wc_xgx,R,'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag);