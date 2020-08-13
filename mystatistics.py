def permutaion_test(df1, df2, function, iter_num = 1000, fraction=0.1):
    buckle = pd.concat([df1, df2])
    size = int(round(min(df1.shape[0]*fraction, df2.shape[0]*fraction), 0))
    res = []
    for i in range(iter_num):
        df = buckle.sample(frac=fraction)
        rdf1 = df.iloc[:size, :]
        rdf2 = df.iloc[size:, :]
        res.append(function(rdf1) - function(rdf2))
    return res
