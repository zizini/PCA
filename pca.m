clear
clc

%We need only masks folder.
listing=dir('C:\Users\User\Desktop\pca\masks\*img'); %Load *.img files.
listing(:).name; %File names.
map=[];

%For each file
for i=1:length(listing) 
    filename=strcat('C:\Users\User\Desktop\pca\masks\',listing(i).name); 
    [avw, machine] = avw_img_read(filename,'','ieee-be'); %Read MRI.
    my_img=avw.img;
    img=my_img(64,:,:); %Select image slice Z=64.
    img=flipud(permute(img,[3,2,1])); %Change orientation to view properly. 
    sz=size(img); %Image's size
    
    mask=(img==512); %We isolate hippocampus. If hippocampus then mask==1 else mask==0.
        
    perimeter=bwperim(mask);%Returns a binary image that contains only the perimeter pixels of hippocampus in the input image mask. Usefull for singed distance map.
     
    signed_distance_map=zeros(sz); %Initialization of singed distance map.
        
    [row col]=find(perimeter); %Find the co-ordinates of the perimeter pixels(one-valued pixels).
      
    %Signed distance map's constructor. 
    %For each signed distance map of size sz.
    %Co-ordinates(x,y)=(rows,colums).
    for x=1:sz(1)
        for y=1:sz(2)
            for j=1:length(row) %For every one-valued pixels.
                temp(j)=sqrt((row(j)-x)^2+(col(j)-y)^2); %Euclidean distance
            end            
            if mask(x,y)%If pixel is inside perimeter set it -, else set it +.  
                signed_distance_map(x,y)=-min(temp);
            else
                signed_distance_map(x,y)=min(temp);
            end
        end
    end   
    
    %After we calculate for each image the signed distance map, it is stored as a column in map array.
    map(:,end+1)=signed_distance_map(:)';
end

[m n]=size(map); 
M=mean(map,2); %For each map's row == for each feature.
imgM=reshape(M,sz); %Reshapes M using the size vector, sz, to define size(imgM). 
perimeterM=bwperim(imgM<0);%Returns a binary image that contains only the perimeter pixels using pixel inside perimeter.
     
M2=map-repmat(M,1,n); %returns an array containing copies of [singed distance map - mean of each feature] with dimensions 1 row and n columns. 
S=M2'*M2/(n-1); %Covariance matrix.

[V D]=eig(S); %returns diagonal matrix D of eigenvalues and matrix V whose columns are the corresponding right eigenvectors of covarience matrix, so that A*V = V*D.
D=diag(D); %The diagonal of D matrix of eigenvalues.

[B,I]=sort(-1*D); %I is the same size as A and describes the arrangement of the elements of -1*D into B along the sorted dimension.
D2=D(I);
V=V(:,I);

D2=cumsum(D2./sum(D2));%Normalization(divide by the sum so as to the sum be 1).
%Cumulative sum: if ë the eigenvalues: ë1, ë2, ë3 etc. comsum => ë1, ë1+ë2, ë1+ë2+ë3. 

D3=D2>=0.97; %Random decision
eigenvalues=find(D3);%Save cumulative sum's of eigenvalues which are >= 0.97. Can be ë1+ë2, ë1+ë2+ë3 ect.
eigenvalues=eigenvalues(1);%Keep only the first

for d=-2:.2:2 %For every sigma.
    for i=1:eigenvalues %For every eigenvalue.    
        eigenvector=V(:,i); %The i-th eigenvectors
        eigenvector=map*eigenvector;%Be part of dataset.
        eigenvector=eigenvector/norm(eigenvector); %Normalization using the 2-norm.
        image=eigenvector*sqrt(D(end-i+1))*d+M; %(eigenvectors)*(index of dispersion d)*(sigma,ë=sigma^2)+(mean image M).
        
        img2=reshape(image,sz); %Devectorize image. 
        new_image=bwperim(img2<0); %Returns a binary image that contains only the perimeter (zero-valued pixels).

        % GIF
        CCall=figure(10); %title(+++); ylabel(+++); xlabel('++++');
        clf
        imshow(new_image,[],'InitialMagnification',500,'Border','tight');%hold on; 
        frame=getframe(CCall);  im=frame2im(frame); [imind,cm]=rgb2ind(im,256);
        filename=strcat('C:\Users\User\Desktop\pca\Output',int2str(i),'.gif');
        if(d==-2) % 1st slide (file creation)
            imwrite(imind,cm,filename,'gif', 'Loopcount',10,'DelayTime',0.1);
        else % rest of gif's slides
            imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.1);
        end
    end
end