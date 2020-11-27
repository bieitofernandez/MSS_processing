function  Sn = FP07noise(params,fr)
    b = params(1);
    m = params(2);
    Sn=(10.^b)*fr.^m;
end