function stop = pswoutfun(optimValues,state)

persistent hist
stop = false;
switch state
    case 'init'
        hist = [0,optimValues.bestfval];
    case 'iter'
        nextline = [optimValues.iteration,optimValues.bestfval];
        disp(nextline)
        hist = [hist;nextline];
    case 'done'
        assignin('base','hist',hist);
end