clear all;
close all;

% Histogram Equalization

Input = imread('D:\University Of Liberal Arts Bangladesh\9th Semester\CSE 429\Assignment\Input.jpg');
J= rgb2gray(Input);

    figure;
        subplot(2,2,1); imshow(J); title('Input Image');
        subplot(2,2,2); imhist(J);title('Input Image Histogram');
        imagehistogram= imadjust(J,[0.3,0.6],[0.0,1.0]);
        newimagehistogram=histeq(J);
        subplot(2,2,3); imshow(newimagehistogram); title('Equalized Image');
        subplot(2,2,4); imhist(newimagehistogram);title('Equalized Image Histogram');
        
        
% Histogram Modeling

K=J;
[M N]=size(K);
K=double(K);

Z=log(1+K);
Zmax=max(Z(:))
Zmin=min(Z(:))

    for i = 1:M
        for j = 1:N
            Y(i,j) = 255*(Z(i,j)-Zmin)/(Zmax-Zmin);
        end
    end
    
    figure();
        subplot(2,2,1);
        imshow(Input);title('Input Image');
        subplot(2,2,2);
        imhist(Input);title('Input Image Histogram');
        subplot(2,2,3);
        imshow(Y,[]); title('Histogram Modeling Image');
        modeling=imread('D:\University Of Liberal Arts Bangladesh\9th Semester\CSE 429\Assignment\modeling.png');
        subplot(2,2,4);imhist(modeling); title('Modeling Histogram');
        
% Histogram Specification
        
a=Input;
a = rgb2gray(a);
b=size(a); 

c=zeros(1,256);

    for i=1:b(1)                                       
        for j=1:b(2)                                    
            for k=0:255                                 
                if a(i,j)==k                            
                    c(k+1)=c(k+1)+1;                    
                end
            end
        end
    end

pdf=(1/(b(1)*b(2)))*c;

cdf = zeros(1,256);
cdf(1)=pdf(1);

    for i=2:256
        cdf(i)=cdf(i-1)+pdf(i);
    end
    
cdf = round(255*cdf);

a1=imread('D:\University Of Liberal Arts Bangladesh\9th Semester\CSE 429\Assignment\srk1.jpg');
a1 = rgb2gray(a1);
b1=size(a1);
a1=double(a1);

c1=zeros(1,256);

    for i1=1:b1(1)
        for j1=1:b1(2)
            for k1=0:255
                if a1(i1,j1)==k1
                    c1(k1+1)=c1(k1+1)+1;
                end
            end
        end
    end
    
pdf1=(1/(b1(1)*b1(2)))*c1;
cdf1 = zeros(1,256);
cdf1(1)=pdf1(1);

    for i1=2:256
        cdf1(i1)=cdf1(i1-1)+pdf1(i1);
    end

cdf1 = round(255*cdf1);

d = 255*ones(1,256);

    for k=1:256
        for k1=1:256
        if cdf(k)<cdf1(k1)
            d(k)=k1;                                   
            break
        end
        end
    end
    
ep = zeros(b(1),b(2));

    for i=1:b(1)
        for j=1:b(2)
            t=(a(i,j)+1);
            ep(i,j)=d(t);
        end
    end


c2 = zeros(1,256);

    for i1=1:b1(1)
        for j1=1:b1(2)
            for k1=0:255
                if ep(i1,j1)==k1
                    c2(k1+1)=c2(k1+1)+1;
                end
            end
        end
    end

    figure();
        subplot(3,2,1);
        imshow(uint8(a));
        title('Input image');
        subplot(3,2,2);
        imhist(Input);
        title('Input image histogram');
        subplot(3,2,3);
        imshow(uint8(a1));
        title('Image with required histogram');
        subplot(3,2,4);
        plot(c1);
        title('Required histogram');
        subplot(3,2,5);
        imshow(uint8(ep));
        title('Modified image');
        subplot(3,2,6);
        plot(c2);
        title('Histogram of modified image');