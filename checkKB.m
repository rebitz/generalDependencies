waiting = 1;

while waiting
    [keyIsDown,secs,keyCode]=KbCheck();
    find(keyCode)
    if keyIsDown==1
        break;
        keyCode(space)
        find(keyCode)
        waiting = 0;
    end

end