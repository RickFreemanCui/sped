def bitsec_lin(m,k,w):
    d = m / w
    cost_mat = (m-k)^3
    num = binomial(d-1, d-1-floor((m-k)/w))
    den = binomial(d-2,d-1-floor((m-k)/w))
    prob_inv = ( num / den )^((1-1/d)*w)
    return log(cost_mat * prob_inv, 2)

def bitsec_isd(m,k,w):
    f = lambda q : q * (1 - log(m/w,2)/(2*m/w)) + log(q / (m-k-q), 2) - k * log(m/w,2) / (2*m/w)
    low, high = 1, m-k-1
    mid = (high-low)/2 + low + 1
    while mid != high:
        if f(mid) > 0:
            high = mid
        elif f(mid) < 0:
            low = mid
        else:
            break
        mid = floor(low + (high - low) / 2 + 1)
    q = mid
    cost1 = q * (k + q) / (m / w) * (m / w)^((k+q)/12)
    cost2 = (m / w)^((k+q)/(m/w)) / 2^q * (m - k - q) * (k+q)/(m/w)
    return log(cost1+cost2, 2)
    
def bitsec_gba(m,k,w):
    r = k - w
    cost = (m / w)^(r/(2*m/w)-1) * r * k
    return log(cost, 2)

def comm_cost (m, k, w):
    return k * (1 - w / m)

def bitsec (m,k,w):
    return min (bitsec_lin(m, k, w), 
                bitsec_isd(m, k, w), 
                bitsec_gba(m, k, w) )
def find_k (params):
    csp = params[0]
    d   = params[1]
    w   = params[2]
    results = params[3]
    m   = w * d
    khigh = floor(m - w * log(d,2))
    klow = ceil (m - log( binomial (m, w), 2 ) )
    min_cost = 999999999999
    min_k = -1
    if bitsec(m,khigh,w) < csp:
        return
    for k in [klow..khigh]:
        if k % d != 0:
            continue
        security = bitsec(m,k,w)
        if security >= csp and comm_cost(m,k,w) < min_cost:
            min_cost = comm_cost(m,k,w)
            min_k = k
            break
    if min_cost != 999999999999:
        results[(csp, d, w)] = min_k

csp_list = [128, 143, 207, 272]
results = dict()
for csp in csp_list:
    for d in [2..10]:
        for w in [100..500]:
            results [(csp, d, w)] = None

for csp in csp_list:
    min_cost = 999999999999
    min_param = (-1, -1, -1)
    for d in [2..10]:
        for w in [100..500]:
            find_k([csp, d, w, results])
            if results[(csp, d, w)] is not None:
                k = results[(csp, d, w)]
                m = d * w
                if comm_cost(m,k,w) < min_cost:
                    min_cost = comm_cost(m,k,w)
                    min_param = (m,k,w)
    m,k,w = min_param
    print (csp, m, k, w, bitsec(m, k, w).n(20), sep='\t')
