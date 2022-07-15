clc;clear; %çalışma alanını temizle

[file,path] = uigetfile("*.jpg","kuşbakışı model seçimi");
secilenResim = fullfile(path,file); % tarama bölgesi dosyasını aç

[file,path] = uigetfile("*.tif","yükseklik model seçimi");
secilenDEM = fullfile(path,file); % tarama bölgesi dosyasını aç


GercekPixelAlani = 0.05*0.05; %m2


resim = imread(secilenResim);%dosyayı oku ve görüntüle

[DEM,DEMR] = readgeoraster(secilenDEM);
for i=1:DEMR.RasterSize(1,1)
    for j=1:DEMR.RasterSize(1,2)
        if DEM(i,j) == -9999
            DEM(i,j) = 0;
        end
        
    end
end

imshow(resim);
hold on


i = 1; endPoint = 0;X = [];Y=[];zm=0;zmctr=0; %nokta seçimi için i değişkeni ve son nokta belirleyici değişken

while i>0
    if endPoint == 0
        [Zx,Zy,button] = myginput(1,"crosshair");
        
    else
        [Zx,Zy,button] = myginput(1,"circle");  
    end
    
    Zx = fix(Zx); Zy=fix(Zy);%resimden alınan x ve y koordinatlarını tam sayı şeklinde değişkene aktarma
 
    
    if(button == 1)
        if endPoint == 1
            mesafeX = abs(Zx - X(1,1));% eğer son nokta işaretçisi açıksa son noktaya olna mesafeyi hesaplama
            mesafeY = abs(Zy - Y(1,1));
        end
        if endPoint==1 && mesafeX < 200 && mesafeY < 200 && i>4 % son noktaya yeteri kadar yakınsak ve tıklama yapıldıysa alanı kapat
            X(i,1) = X(1,1);
            Y(i,1) = Y(1,1);
            i=-1;
            fill(X,Y,"g","FaceAlpha","0.3","LineStyle","none");% alan seçimi bitince alanı doldur
        else%değilse nokta eklemeye devam et
            X(i,1) = Zx;
            Y(i,1) = Zy;
        end 
        i=i+1;
    end
    %%hatalı nokta silme
    
    if button == 127 && ~isempty(X) %% son nokta DEL tuşuna basılırsa silinecek
        tempx = [];tempy = [];
        for j=1:1:i-2
            tempx(j,1) = X(j,1);
            tempy(j,1) = Y(j,1);
        end
        X = [];
        Y= [];
        X = tempx;Y=tempy;
        clf("reset");
        imshow(resim);
        hold on;
        i= i-1;
    end
    
    %%uzaklaştırma yaklaştırma
    if button == 43 % + butonuna basılınca yakınlaştır 
        zmctr = zmctr +1;
        if zm == 0
            sol = Zx - 300;
            sag = Zx + 300;
            ust = Zy - 300;
            alt = Zy + 300;
            zm=1;
        elseif ~(abs(sol-sag)<101 && abs(ust-alt)<101)% yakınlaştırma oranını aşmamak için kontrol et
            sol = sol + 50;
            sag = sag - 50;
            ust = ust + 50;
            alt = alt - 50;
        end
        xlim([sol,sag]);ylim([ust,alt]);
        
    end
    if button == 45
       if zmctr > 0 && zmctr < 3
           sol = 1;
           sag = DEMR.RasterSize(1,2)-1;
           ust = 1;
           alt = DEMR.RasterSize(1,1)-1;
           zmctr = 0;
           zm = 0;
       elseif zmctr>=3
           sol = sol - 50;
           sag = sag + 50;
           ust = ust - 50;
           alt = alt + 50;
           zmctr  = zmctr - 1;
           zm=0;
            
       end
       xlim([sol,sag]);ylim([ust,alt]);
    end
     
        
    
    %işaretlemeler için görsellik
        plot(X,Y,"--bx",...
        "LineWidth",2,...
        "MarkerSize",10,...
        "MarkerEdgeColor","r",...
        "MarkerFaceColor",[0.5,0.5,0.5]);
    
    
    
    % S tuşuna basıldığında son nokta bulucuyu aktifleşitrme/pasifleştirme
    if button == 115 && endPoint == 0
        endPoint = 1;
    elseif button == 115 && endPoint == 1
        endPoint = 0;
    end
    
    
    %programı ESC tuşu ile sonlandır
    if button == 27
        clf;
        i=-1;
    end
end

%inpolygon fonksiyonunu çağırabilmek için satır ve sütun numaralarını kedi
%kendine atama yapıyoruz
Koord(:,:,1) = zeros(DEMR.RasterSize(1,1),DEMR.RasterSize(1,2));
Koord(:,:,2) = zeros(DEMR.RasterSize(1,1),DEMR.RasterSize(1,2));


msg2 = waitbar(0,"KOORDİNAT VERİLERİ YERLEŞTİRİLİYOR");
msg1 = waitbar(0,"VERİLER ALINIYOR");

for i=1:DEMR.RasterSize(1,1)
    waitbar(i/DEMR.RasterSize(1,1),msg1);
    Koord(i,:,2) = i;
end
close(msg1);
for i=1:DEMR.RasterSize(1,2)
    waitbar(i/DEMR.RasterSize(1,2),msg2);
    Koord(:,i,1) = i;
end
close(msg2);

[in,on] = inpolygon(Koord(:,:,1),Koord(:,:,2),X,Y);

temp = 1,realValues = [];
for i=1:DEMR.RasterSize(1,1)
    for j=1:DEMR.RasterSize(1,2)
        if in(i,j) == 1
           Yson(temp,1) = i;
           Xson(temp,1) = j;
            
           realValues(temp,1) = DEM(Yson(temp,1),Xson(temp,1));
           temp = temp +1;
        end
    end
end

zemin = min(realValues); % zemin bilgisi girildi

toplamHacim = 0;
for i=1:length(realValues)
    farklar(i,1) = realValues(i,1) - zemin;%yükseklik farklarını zemin referansına göre tek tek hesaplama işlemi
    toplamHacim = toplamHacim + GercekPixelAlani * farklar(i,1);
end


realValuesForGraph = [];
for i=1:length(farklar)
    realValuesForGraph(i) = farklar(i)*50;
end
    imshow(resim);%resmin üzerinde 3d gösterme
    hold on;
  
    xlabel("TOPLAM HACİM: " + toplamHacim + "m3")
    plot3(Xson,Yson,realValuesForGraph,"o",... 
        "MarkerSize",0.5,...
        "MarkerEdgeColor","g",...
        "MarkerFaceColor",[0.5,0.5,0.5]);




