function stop=psooutfcn(optimValues,state)
persistent X
stop=false;
switch state
case 'init'
        X= optimValues.swarm;
    case 'iter'
    nextline = optimValues.swarm;
        X = [X;nextline];
    case 'done'
        assignin('base','X',X);

end