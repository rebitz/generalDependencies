letters = {'a','b','c','d','e','f','g','h','i','j','k','l','m','n','o',...
    'p','q','r','s','t','u','v','w','x','y','z'};

nLetters = 3; % letters in the code
nToGen = 50;


tmp = letters(randi(length(letters),nToGen,nLetters));

fu = [];
for i = 1:nToGen
    fu = [fu; strcat(tmp{i,:})]
end