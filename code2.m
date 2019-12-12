clc
clear
%%
% M=190
Add=double(zeros(31266,1));
% For loop to go through 190 images and vetorize them. In the same loop I
% have used the Add variable to add all the images. The 'Add' variable will
% be used later to find the average of all the images and stored in
% variable 'psi'.
for i=1:190  
    im=double(imread(strcat(int2str(i),'a.jpg')));
    im=reshape(im',[31266,1]);
    gamma(:,i)=im;
    Add=Add+im;
end
%% Average Image
psi=(Add/190);
psi_img=uint8(reshape(psi,[162,193]))'; % Reshaping 'psi' to 193x162 to 
imshow(psi_img);                        % display  Average image 
title('Average Image')
%% Getting the A matrix by subtrating the Average face vector 'psi' from the 
%% each 190 vectorized image and storing them in the columns of A.
for i=1:190
   phi(:,i)=(gamma(:,i))-psi;
end
A=phi;
A=double(A);
%% Plot singular values of Matrix B
L=((A')*A); %  Calculating the L matrix i.e. L = (A^T)*A
[U,S,V]=svd(A);
% V is the eigenvector matrix of A'A or L
% Finding the Eigenvalues and Eigenvectors of L using MATLAB function eigs, 
% here V- Eigen Vector matrix, D- Diagonal Matrix of Eigen Values sorted 
% in descending order, i.e. D(1,1)>D(2,2)>D(3,3) ...... D(190,190) 
%[V,D]=eigs(L,190);
 % [V,D]=eigs(B,70);

% Singular Value Plot for the Eigenvalues of data matrix

for i=1:190
    sing_Val(i)=S(i,i); % Putting Diagonal Values of S into eigvalue array. 
end
figure,plot((sing_Val)) % Plotting the singular values of the data in descending order
title('Singular Value plot of the Data');
ylabel('Singular Value magnitude');
%% Section to find the most important PCs.
eigvalue=sing_Val.^2; % Getting the eigenvalues from singular values.
eigsum=sum(eigvalue);
csum=0;
k90=0;
% d=10;
% for d=1:190
    for i=1:190
       csum=csum+(eigvalue(i));
       tv=csum/eigsum;
       if tv>0.95
         k90=i;
         break;
       end
    end
% end
fprintf('No of Principal components: %d \n',k90)

%% Getting the Eigen vector matrix for (A*(A^T)) matrix. 
%% In other words finding the eigen Face vectors
eignVect_A=A*V; % 

%% Reconstructing 107a.jpg image using different nos of PCs. PC=1:190
%% finally plotting the MSE vs no of PCs plot
mse_normal=[];
imge=double(imread(strcat('107a.jpg')));
imge=reshape(imge',[31266,1]);
w1=[];
for n=1:190 % For loop to plot the MSE vs No of PCs variation. Here n= No of PCs
    Vmax=V(:,1:n); %getting the correspoding "orthonormal Eigenvectors" or PCs
    eignVect_phi_Max=A*Vmax;
    % Now I normalized the calculated Eigenvectors, eslse some pixel values
    % exceed 255 and as a result proper reconstruction of the image doesn't
    % occur.
    for i=1:n
        eignVect_phi_Max(:,i)=eignVect_phi_Max(:,i)/norm(eignVect_phi_Max(:,i));
    end
    % Calculating the weight vector w where each element of w represents
    % the contribution of each eigenface vector to the reconstructed image.
    w1= (eignVect_phi_Max')*(imge-double(psi)); 
    rec_img=(eignVect_phi_Max*w1); % Reconstructing the image
    err=imge-(rec_img+psi); % Calculating the error w.r.t original image
    mse_normal(n)=(((err')*err)/31266); % Calculting Mean Square Error.
    
    % Displaying the reconstruction process as the number of PCs increases.
    % Uncomment the line below to see the reconstruction process.
    %imshow((reshape(uint8(rec_img),[162,193]))'+psi_img);   // 
end

% %Plotting the final reconstructed image and MSE w.r.t no of PCs.
figure, plot(mse_normal);
title('MSE varaiation plot for Normal Face as No of PCs are varied');
xlabel('No. of Principal Components');
ylabel('Mean Square Error');
figure,subplot(1,2,1);
imshow(imread(strcat('107a.jpg')));
title('Original Image');
subplot(1,2,2);
% figure
imshow((reshape(uint8(rec_img),[162,193]))'+psi_img);
title('Reconstructed Image with 190 PCs');

%% MSE for Smiling images

imge=double(imread(strcat('107b.jpg')));
imge=reshape(imge',[31266,1]);
w1=[];
% figure;
for n=1:190
% n=70;
Vmax=V(:,1:n); %getting the correspoding "orthonormal Eigenvectors" or PCs
    eignVect_phi_Max=A*Vmax;
    % Now I normalized the calculated Eigenvectors, else some pixel values
    % exceed 255 and as a result proper reconstruction of the image doesn't
    % occur and we get white patches here and there on the image.    
    for i=1:n
        eignVect_phi_Max(:,i)=eignVect_phi_Max(:,i)/norm(eignVect_phi_Max(:,i));
    end
    % Calculating the weight vector w where each element of w represents
    % the contribution of each eigenface vector to the reconstructed image.
    w1= (eignVect_phi_Max')*(imge-psi);
    rec_img=(eignVect_phi_Max*w1);% Reconstructing the image
    err=imge-(rec_img+psi);% Calculating the error w.r.t original image
    mse_smile(n)=(((err')*err)/31266);% Calculting Mean Square Error.
    % Displaying the reconstruction process as the number of PCs increases.
    % Uncomment the line below to see the reconstruction process.
    %imshow((reshape(uint8(rec_img),[162,193]))'+psi_img);
end
figure, plot(mse_smile);
% title('MSE Smiling Face')
title('MSE varaiation plot for Smiling Face as No of PCs are varied');
xlabel('No. of Principal Components');
ylabel('Mean Square Error');
figure,subplot(1,2,1);
imshow(imread(strcat('107b.jpg')));
title('Original Image');
subplot(1,2,2);
imshow((reshape(uint8(rec_img),[162,193]))'+psi_img);
title('Reconstructed Image with 190 PCs');

%% MSE for normal image in the test set

imge=double(imread(strcat('192a.jpg')));
imge=reshape(imge',[31266,1]);
w1=[];
% figure;

for n=1:190
% n=70;
Vmax=V(:,1:n); %getting the correspoding "orthonormal Eigenvectors" or PCs
    eignVect_phi_Max=A*Vmax;
    % Now I normalized the calculated Eigenvectors, else some pixel values
    % exceed 255 and as a result proper reconstruction of the image doesn't
    % occur.  
    for i=1:n
        eignVect_phi_Max(:,i)=eignVect_phi_Max(:,i)/norm(eignVect_phi_Max(:,i));
    end
    % Calculating the weight vector w where each element of w represents
    % the contribution of each eigenface vector to the reconstructed image.
    w1= (eignVect_phi_Max')*(imge-psi);
    rec_img=(eignVect_phi_Max*w1);% Reconstructing the image
    err=imge-(rec_img+psi);% Calculating the error w.r.t original image
    mse_norm_testset(n)=(((err')*err)/31266);% Calculting Mean Square Error.
   % Displaying the reconstruction process as the number of PCs increases.
    % Uncomment the line below to see the reconstruction process. 
%     imshow((reshape(uint8(rec_img),[162,193]))'+psi_img);
%     title('Reconstructed image')
end
% subplot(1,2,1);
figure, plot(mse_norm_testset);
title('MSE varaiation plotfor Normal face from Test Set');
xlabel('No. of Principal Components');
ylabel('Mean Square Error');
figure,subplot(1,2,1);
imshow(imread(strcat('192a.jpg')));
title('Original Image');
subplot(1,2,2);
imshow((reshape(uint8(rec_img),[162,193]))'+psi_img);
title('Reconstructed Image with 190 PCs');
%% Reconstructing Non-Human image (image of flower Vase) with all the eigenvectors/PCs

imge=double(imread('vase1.jpg'));
imge=double(reshape(imge',[31266,1]));
% The next loop is to make sure the eigenvectors are normalized. 
% Although not necessary, just a precaution
for i=1:190 
        eignVect_A(:,i)=eignVect_A(:,i)/norm(eignVect_A(:,i));
end
    % Calculating the weight vector w where each element of w represents
    % the contribution of each eigenface vector to the reconstructed image.
w_test=(eignVect_A')*(imge-psi); % Weights

% Reconstruction
rec_img=uint8(eignVect_A*w_test); % Reconstructing the image
q=(reshape(uint8(rec_img),[162,193]))'+psi_img; % Resizing the image

% Displaying the images
figure,subplot(1,2,1);
imshow(imread(strcat('vase1.jpg')));
title('Original Image');
subplot(1,2,2);
imshow(q);
title('Reconstructed Image with 190 PCs');
%% Reconstructing Rotated Image by rotating it to different angles and using all the PCs. 
%% Here I rotated the image by 1 degree anticlockwise everytime and plotted the MSE curve.
% I have also used this code snippet to also manually put values of angle of
% rotation and get the reconstructed image.

% The next loop is to make sure the eigenvectors are normalized. 
% Although not necessary, just a precaution
for i=1:190
        eignVect_A(:,i)=eignVect_A(:,i)/norm(eignVect_A(:,i));
end
%%
%  for k=1:360 %% uncomment this line, line 222,230,231,233,235and 238 to 241.
          % to get MSE plot for every 1 degree rotation of the input image.
%      J=imrotate(imread('120a.jpg'),k,'bilinear','crop');
     J=imrotate(imread('120a.jpg'),355,'bilinear','crop'); % comment this to get mse plot
     imge=double(reshape(J',[31266,1]));    
    
     w_test=(eignVect_A')*(imge-psi); % Weights

% Reconstruction
     rec_img=((eignVect_A*w_test)); % Reconstructing the image
%     err=imge-(rec_img+psi); % Calculating the error w.r.t original image
%     mse_rotate(k+1)=(((err')*err)/31266); % Calculting Mean Square Error.
     q=(reshape(uint8(rec_img),[162,193]))'+psi_img; % Resize the image
%     imshow(q); % show the live plot.
     
%  end
 %%
% figure,plot(mse_rotate)
% title('MSE Plot for rotated imaage (0-360)');
% xlabel('Input image rotation angle ( in Degrees) ')
% ylabel('Mean Squared Error');
 %%
figure,subplot(1,3,1);
imshow(imread(strcat('120a.jpg')));
title('Original Image');
subplot(1,3,2);
imshow(J);
title('Rotated by 355 degree anticlockwise');
subplot(1,3,3);
imshow(q);
title('Final Rotated image reconstruction');
