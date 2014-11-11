function lambda = expDecayLambda(n)
    global lambdaIni lambdaSS gamma nIni
    lambda = lambdaIni / (n+nIni) ^ gamma + lambdaSS;
end