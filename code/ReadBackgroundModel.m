function [B0, B1]=ReadBackgroundModel(backgroundModelFile)

[a,c,g,t]=textread(backgroundModelFile,  ...
                   '%f %f %f %f', ...
                   'commentstyle', 'shell', ...
                   'emptyvalue', nan);
                   

B0 = [a(1) c(1) g(1) t(1)];
B1 = [a(6:9) c(6:9) g(6:9) t(6:9)];


                   
                  