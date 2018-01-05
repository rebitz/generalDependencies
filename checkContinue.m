% check for continuation, useful for human subjects, to allow them to
% respond with a button press to go beyond the next stage of information

function next = checkContinue(key)

    curTime = GetSecs;
    
    next = 0;
    
    [keyIsDown,secs,keyCode] = KbCheck;
            
    % Key pressed -> fixation acquired
    if keyIsDown && keyCode(key) && secs > curTime;
        next = 1;
    else
        next = 0;
    end
    
    simple_esc_check;
    
end
    