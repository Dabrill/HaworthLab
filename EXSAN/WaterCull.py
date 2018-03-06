recorder = None
def PRINT(string):
    global recorder
    print(string)
    if recorder is not None:
        recorder.write(string+"\n")
def setRecorder(obj):
    global recorder
    recorder = obj
def stochasticWaterCull(CONSTANTS,step,data):
    limTemp = CONSTANTS["limTemp"]
    limTemp = step.get("limTemp",limTemp)
    limPoses = CONSTANTS["limPoses"]
    limPoses = step.get("limPoses",limPoses)
    limExp = CONSTANTS["limExp"]
    limExp = step.get("limExp",limExp)
    padding = CONSTANTS["padding"]
    padding = step.get("padding",padding)
    TEMP = limTemp * limPoses**limExp / data.numGood()**limExp
    PRINT("Calculated \"Temperature\" is "+str(TEMP))
    PRINT("Current Pose Count "+str(data.numGood())+"/"+str(data.getLength()))
    poses = data.fixvar.poses
    best = data.fixvar.bestPose
    bestScore = best.energy
    PRINT("BEST: "+str(bestScore))
    kept = 0
    outf = open("StochProbTable.txt",'w')
    for p in poses:
        if p.energy is not None:
            roll = random.random()
            chance = math.exp((padding*bestScore)/TEMP) / math.exp(p.energy/TEMP)
            #PRINT(str(p.energy)+" - "+str(chance))
            outf.write(str(p.energy)+"\t"+str(chance)+"\t"+str(p.retain)+"\n")
            if (roll > chance):
                data.bad(p.num)
            else:
                if (p.retain):
                    kept+=1
    outf.close()
    PRINT("\tPoses retained "+str(kept))
    data.saveCullResults("WaterCull.pos") 
def rubberTopWaters(step,data):
    EXPONENT_CONSTANT = 0.35
    limNum = step["liminalNumber"]
    numGood = data.numGood()
    numKeep = min(1,limNum**EXPONENT_CONSTANT/numGood**EXPONENT_CONSTANT)*numGood
    topWaters(step,data,numKeep)
    poses = data.fixvar.poses
def topWaters(step,data,explicitTop = None):
    def veryHighIfNone(x):
        if (x.energy is not None):
            return x.energy
        return 999999
    if explicitTop is None:
        numKeep = min(step["liminalNumber"],data.numGood())
    else:
        numKeep = min(explicitTop,data.numGood())
        limNum = explicitTop
    numGood = data.numGood()
    numTotal = data.getLength()
    PRINT("Current Pose Count "+str(data.numGood())+"/"+str(data.getLength())+" of which "+str(numKeep)+" will be kept")
    poses = [None] * numTotal
    for p in range(numTotal):
        poses[p] = data.fixvar.poses[p]
    outf = open('topwaters.dtf','w')
    poses = sorted(poses,key=veryHighIfNone)
   
    num = 0 
    cur = 0
    while (cur < numKeep):
        if (poses[num].retain):
            outf.write(str(poses[num].num)+"\n")
            cur+=1
        num+=1
    while (num < numTotal):
        data.bad(poses[num].num)
        num+=1
    outf.close()
    data.saveCullResults("WaterCull.pos")
def thresholdOfMaxEnergy(step,data):
    frac = step["threshold"]
    isRelative = step.get("relative",False)
    PRINT("Creating table of poses with at least %3.1f%% of max energy"%(frac*100))
    PRINT("Current Pose Count "+str(data.numGood())+"/"+str(data.getLength()))
    
    num = 0
    poses = data.fixvar.poses
    best = data.fixvar.bestPose
    worst = best
    if isRelative:
        floorPose = best
        for p in poses:
            if (p.retain) and (p.energy is not None):
                if (p.energy > floorPose.energy):
                    floorPose = p
        PRINT("Floor has energy "+str(floorPose.energy)+" from pose "+str(floorPose.num))
        floor = floorPose.energy
    else:
        floor = 0
        PRINT("Floor energy is default zero (0)")
    liminal = ((best.energy - floor) * frac) + floor
    PRINT("Liminal energy is "+str(liminal))
    for p in poses:
        if (p.energy is not None):
            if (p.energy > liminal):
                data.bad(p.num)
            else:
                num+=1
                if p.retain and (p.energy > worst.energy):
                    worst = p
    PRINT("Best energy "+str(best.energy)+" from pose "+str(best.num))
    PRINT("Lowest included energy "+str(worst.energy)+" from pose "+str(worst.num))
    PRINT("\tPoses retained "+str(num))
    data.saveCullResults("WaterCull.pos")
def getDeltaFrac(topDir,plan,i):
    step = plan[i]
    ceiling = step["threshold"]
    scaling = step.get("scaling",-900)
    prior = None
    now = None
    mode = 0
    while (i >= 0):
        step = plan[i]
        if (step["type"] == "water"):
            mode = 1
        if (step["type"] == "TMD"):
            if (mode == 1):
                if (now is None):
                    now = step["folder"]
                elif (prior is None):
                    prior = step["folder"]
            mode = 0
        i-=1
    delta = averageEnergyChange("%s/%s"%(topDir,prior),"%s/%s"%(topDir,now))
    frac = min(ceiling,delta/scaling)
    return frac
def cutEnd(aList):
    endpoint = len(aList) - 1
    while (len(aList[endpoint]) < 1):
        endpoint-=1
    return aList[:endpoint+1]
def makeAlignment(prior,now):
    prior = cutEnd(prior)
    now = cutEnd(now)
    alignment = [None] * len(now)
    j = 0
    nowPose = now[0]
    for i in range(len(prior)):
        priorPose = prior[i] 
        while(priorPose == nowPose[0:len(priorPose)]):
            '''
            print(i,j)
            print("\t%s"%priorPose)
            print("\t%s"%nowPose)
            input("")
            '''

            alignment[j] = i
            j+=1
            if (j < len(now)):
                nowPose = now[j]
            else:
                return alignment
    return alignment
def averageEnergyChange(preFolder,nowFolder):
    from ETable import EnergyTable
    #print("Start")
    prior = open("%s/fixvar.out"%preFolder,'r').read().split("\n")
    prior = prior[3:]
    now = open("%s/fixvar.out"%nowFolder,'r').read().split("\n")
    now = now[3:]
    alignment = makeAlignment(prior,now)
    #print("Fixvar done")
    priorEnergy = EnergyTable.readTableFile("%s/ETable.tsv"%preFolder)
    nowEnergy = EnergyTable.readTableFile("%s/ETable.tsv"%nowFolder)
    #print("Energy done")
    deltaDict = {}
    for i in range(len(alignment)):
        a = alignment[i]
        nowEnergyPose = nowEnergy.getParameter("Energy",i+1)
        
        if nowEnergyPose is not None:
            if a not in deltaDict:
                deltaDict[a] = []
            priorEnergyPose = priorEnergy.getParameter("Energy",a+1)
            delta = nowEnergyPose - priorEnergyPose
            #print(a,priorEnergyPose,nowEnergyPose,delta)
            #input("")
            deltaDict[a].append(delta)
    #print("Grouped")
    grandTotal = 0.0
    grandCount = 0
    for key in deltaDict:
        entry = deltaDict[key]
        #total = sum(entry)
        #mean = total / len(entry)
        #entry = sorted(entry)
        #median = entry[len(entry)//2]
        highest = min(entry)
        #print("Alignment for prior pose %i Members: %i Mean: %4.1f Median: %4.1f Max: %4.1f"%(key,len(entry),mean,median,highest))
        #print(key,"Average ",mean," Median ",median," highest ",highest)
        grandTotal+=highest
        grandCount+=1
    return grandTotal/grandCount        
    
