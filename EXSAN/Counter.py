import time
import sys
class Counter:
    def __init__(me,n,template = "Poses: %i/%i"):
        me.last = time.time()
        me.n = n
        me.idleMode = ("idlelib" in sys.modules)
        me.template = template
    def disp(me,i):
        now = time.time() 
        if me.idleMode:
            '''
            if (i % 1000 == 0):
                print("%i/%i"%(i,me.n))
            '''
            if (now > me.last + 10):
                me.last = now
                print(me.template%(i,me.n))
        else:
            #now = time.time() 
            if (now > me.last + 0.2):
                me.last = now
                print(me.template%(i,me.n),end="\r")
            if (i+1 == me.n):
                print(me.template%(me.n,me.n),end="\r")
                print("\nDone!")
def fakeFunction():
    n = 140000
    st = time.time()
    cnt = Counter(n)
    for i in range(n):
        cnt.disp(i)
        if (i % 29 == 0):
            time.sleep(.001)
    print((time.time()-st))
#fakeFunction()
#time.sleep(100)
