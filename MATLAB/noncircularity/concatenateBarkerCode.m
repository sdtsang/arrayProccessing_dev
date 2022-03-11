% Calculate concateated Barker Codes Correctly:

% To generate a barker code with length 15, "cascade" barker code in which
% you can build longer barker code sets from shorter sets

% 1.) [1 1 1 -1 1] x [1 1 -1] = [[1 1 1 -1 1][1 1 1 -1 1][-1 -1 -1 1 -1]] 
% barker_3 * barker_5 = barker_15
% barker_15 = [1 1 1 -1 1 1 1 1 -1 1 -1 -1 -1 1 -1];

% 2.) [1 1 1 -1 1] x [1 1 1 -1 -1 1 -1] = barker_5 x barker_7 = barker_35

% 3.) [1 1 1 -1 -1 1 -1] x [1 1 1 -1 -1 -1 1 -1 -1 1 -1] = barker_7 x
% barker_11 = barker_77

% 4.) [+1 +1 +1 −1 −1 −1 +1 −1 −1 +1 −1] x [+1 +1 +1 +1 +1 −1 −1 +1 +1 −1
% +1 −1 +1] = barker_11 x barker_13 = barker_143

% write for loops to generate each:
barker_3 = [1 1 -1];
barker_5 = [1 1 1 -1 1];
barker_7 = [1 1 1 -1 -1 1 -1];
barker_11 = [1 1 1 -1 -1 -1 1 -1 -1 1 -1];
barker_13 = [+1 +1 +1 +1 +1 -1 -1 +1 +1 -1 +1 -1 +1];

% test and verify on barker 15:
barker_15 = zeros(length(barker_5),length(barker_3));
for i = 1:length(barker_3)
    barker_15(:,i) = barker_5 * barker_3(i);
end
barker_15 = barker_15(:);   

% barker 35
barker_35 = zeros(length(barker_7),length(barker_5));
for i = 1:length(barker_5)
    barker_35(:,i) = barker_7 * barker_5(i);
end
barker_35 = barker_35(:);    

% barker_77
barker_77 = zeros(length(barker_11),length(barker_7));
for i = 1:length(barker_7)
    barker_77(:,i) = barker_11 * barker_7(i);
end
barker_77 = barker_77(:);    
    
% barker_143
barker_143 = zeros(length(barker_13),length(barker_11));
for i = 1:length(barker_11)
    barker_143(:,i) = barker_13 * barker_11(i);
end
barker_143 = barker_143(:);
    