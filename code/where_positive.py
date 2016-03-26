def where_positive(x):
    l = len(x)
    i = 0
    intervals = []
    a = None
    while i < l:
        if a is None and x[i] > 0:
            a = i
        if a is not None and (x[i] <= 0 or i == l - 1):
            intervals.append([a,i])
            a = None
        i += 1
    return intervals

# print(where_positive([0,0,0,0,1,1,2,3,4,1,0,1,0,2,1,1,1]))