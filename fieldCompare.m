function fieldCompare(struct1,struct2)

    s1 = fieldnames(struct1);
    s2 = fieldnames(struct2);
    
    fu = ismember(s1,s2);
    if sum(fu) ~= length(s1)
        disp('first struct includes additional field(s):')
        fields = s1(~fu)
    else
        disp('first struct is fully contained within the second.')
    end
    
    fu = ismember(s2,s1);
    if sum(fu) ~= length(s2)
        disp('second struct includes additional field(s):')
        fields = s2(~fu)
    else
        disp('second struct is fully contained within the first.')
    end