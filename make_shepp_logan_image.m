%Shepp Logan Phantom, Fourier and Filtered Rec, Correct BD and Geg Rec:  No Noise

function phantomData = make_shepp_logan_image(imageResolution)

phantomData=zeros(imageResolution);
fourier_n = imageResolution/2;

for j=1:2*fourier_n
   for k=1:2*fourier_n
      x=(j-fourier_n-1)/fourier_n;
      y=(k-fourier_n-1)/fourier_n;
      xi=(x-.22)*cos(.4*pi)+y*sin(.4*pi);
      eta=y*cos(.4*pi)-(x-.22)*sin(.4*pi);
      
      z=0;
      if (x/.69)^2+(y/.92)^2<=1
         z=2;
      end
      if (x/.6624)^2+((y+.0184)/.874)^2<=1
         z=z-.98;
      end
      if (xi/.31)^2+(eta/.11)^2<=1
         z=z-.8;
      end
      
      xi=(x+.22)*cos(.6*pi)+y*sin(.6*pi);
      eta=y*cos(.6*pi)-(x+.22)*sin(.6*pi);
      
      if (xi/.41)^2 +(eta/.16)^2<=1
         z=z-.8;
      end
      if (x/.21)^2+((y-.35)/.25)^2<=1
         z=z+.4;
      end
      if (x/.046)^2+((y-.1)/.046)^2<=1
         z=z+.4;
      end
      if (x/.046)^2+((y+.1)/.046)^2<=1
         z=z+.4;
      end
      if ((x+.08)/.046)^2+((y+.605)/.023)^2<=1
         z=z+.4;
      end
      if (x/.023)^2+((y+.605)/.023)^2<=1
         z=z+.4;
      end
      if ((x-.06)/.023)^2+((y+.605)/.046)^2<=1
         z=z+.4;
      end
      phantomData(j,k)=z;
   end
end
phantomData = transpose(phantomData);
phantomData = phantomData(end:-1:1,:);
