function y=svm2mat(s)
% Usage: y=svm2mat ('filename')
% clc
fid = fopen(s);
i=1;
while ~feof(fid) % not end of the file 
       s = fgetl(fid); % get a line 
       s1=[];  j=1;
       while (j<=length(s))
             while(j<=length(s) && s(j)~=' ')
                   s1=[s1 s(j)];
                   j=j+1;   
             end  
             j=j+1;   
             s1=[s1 ' '];              
             while ( (j<length(s)) && (s(j)~=':') )   
                    j=j+1;   
             end             
             j=j+1;   
       end
        s2=str2num(s1) ;
        if (i==1)
            yy=zeros(1,length(s2));
        end
        yy=[yy ; s2];
        i=i+1;      
end
yy(1,:)=[];
y=yy;
