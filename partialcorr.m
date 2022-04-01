[a,R]=geotiffread('E:\test\npp_ccy\occipNPP_2003_01_01.tif');%先导入投影信息，某个影像的路径就行（最好是你分析的数据中的一个）
info=geotiffinfo('E:\test\npp_ccy\occipNPP_2003_01_01.tif');%同上
[m,n]=size(a);
warning('off');
nppsum=zeros(m*n,18);%此处要修改，共几年就填写多少，我这里是8年的

for year=2003:2020
    filename=strcat('E:\test\npp_ccy\occipNPP_',int2str(year),'_01_01.tif');%此处要修改，我这里是八年的每年年均植被覆盖度的数据，注意你的文件名字。
    data=importdata(filename);
    data=reshape(data,m*n,1);
    nppsum(:,year-2002)=data;%此处需要修改，我的数据是从2013开始，此处就为2012.
end

wcsum=zeros(m*n,18);
for year=2003:2020
    filename=strcat('E:\test\pre_ccy\ocuqi_pre_',int2str(year),'_1.tif');%此处要修改，我这里是八年的每年年均地表温度的数据，注意你的文件名字。
    data=importdata(filename);
    data(m+1,:)=[];
    data(:,n+1)=[];
    data=reshape(data,m*n,1);
    wcsum(:,year-2002)=data;
end


for year=2003:2020
    filename=strcat('E:\test\wendu_units\ocuqi_temp_',int2str(year),'_1.tif');%此处要修改，我这里是八年的每年年均地表温度的数据，注意你的文件名字。
    data=importdata(filename);
    data(m+1,:)=[];
    data(:,n+1)=[];
    data=reshape(data,m*n,1);
    wcsumz(:,year-2002)=data;
end




count = 2020-2003+1;
%相关性和显著性
npp_wc_xgx=zeros(m,n);
npp_wc_p=zeros(m,n);
for i=1:length(nppsum)
    npp=nppsum(i,:);
    if min(npp)>=0 %注意这里的NPP的有效范围是大于0，如果自己的数据有效范围有小于0的话，则可以不用加这个
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

filename5='E:\test\jiangyu_units\ocu\shiyan_pre5p.tif';%此处要修改，输出的路径及名字
filename6='E:\test\jiangyu_units\ocu\显著性_pre25.tif';%同上


%%输出图像
geotiffwrite(filename5,npp_wc_xgx,R,'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag);
geotiffwrite(filename6,npp_wc_p,R,'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag);



npp_wc_xgx(npp_wc_p<2.10092)=NaN;
name1='E:\test\jiangyu_units\ocu\通过显著性0.05检验的年偏相关03-20趋势值05.tif';
geotiffwrite(name1,npp_wc_xgx,R,'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag);