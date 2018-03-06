import math
def gamma(n):
    i = int(n)
    diff = n - i
    if (diff > 0.49):
        m = int(2 * n)
        top = doubleFactorial(m-2)
        bottom = 2**((m-1)/2)
        return top / bottom * math.sqrt(math.pi)
    return factorial(i-1)
def factorial(n):
    a = 1
    for m in range(1,n+1):
        a = m * a
    return a
def doubleFactorial(n):
    at = n
    ans = 1
    while(at > 0):
        ans = ans * at
        at-=2
    return ans
def beta(x,y):
    a = math.sqrt(2*math.pi)
    top = x**(x-0.5)*y**(y-0.5)
    bottom = (x+y)**(x+y-0.5)
    return a * top / bottom
def probDense(t,v):
    a = (v+1)/2.0
    b = v /2.0
    tSquared = (t**2)
    tOverV = tSquared/ v

    left = gamma(a) / (math.sqrt(v*math.pi)*gamma(b))
    right = (1 + tOverV)**(-a)
    return left * right
def studentT(c,dof,tails = 2):
    if (tails < 1) or (tails > 2):
        raise Exception("Illegal number of tails")
    x = -c
    grain = 0.001
    cum = 0
    while (x < c):
        cum+=probDense(x,dof) * grain
        x+=grain
    return round((1 - cum) * (tails / 2),3)
def ttest(L1,L2,tails = 2):
    if (tails < 1) or (tails > 2):
        raise Exception("Illegal number of tails")
    if (len(L1) == 0) or (len(L2) == 0):
        return None
    num1=sum(L1)
    den1=len(L1)
    avg1=num1/den1

    num2=sum(L2)
    den2=len(L2)
    avg2=num2/den2
    if (den1+den2 < 3):
        return None

    sumdiff=0
    for i in L1:
        diff1=((i-avg1)**2)
        sumdiff=sumdiff+diff1

    sumdiff2=0
    for j in L2:
        diff2=(j-avg2)**2
        sumdiff2=sumdiff2+diff2

    varden=den2+den1-2
    varience=(sumdiff+sumdiff2)/varden
    y=math.sqrt((varience*((1/den1)+(1/den2))))
    t=((abs(avg1-avg2))/y)
    df=varden
    return studentT(t,df,tails)
