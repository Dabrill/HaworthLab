import multiprocessing
import Counter
recorder = None
def PRINT(string):
    global recorder
    print(string)
    if recorder is not None:
        recorder.write(string+"\n")
def setRecorder(obj):
    global recorder
    recorder = obj
def distCullThread(inQueue,outQueue,tree,maxCut,resNum):
    while(True):
        next_task = inQueue.get()
        if next_task is None:
            inQueue.task_done()
            break
        tup = (next_task.indexNum,next_task.resInDistance(tree,maxCut,resNum))
        outQueue.put(tup)
        inQueue.task_done()
def distCull(tree,step,poseData,resNum, cores = 1):
    maxCut = step["max"]
    PRINT("Performing distance cull on "+str(poseData.getLength())+" poses")

    poseList,includeAll = poseData.goodPoseList()
    numberOfPoses = len(poseList)
    jobs = []
    tasks = multiprocessing.JoinableQueue()
    results = multiprocessing.Queue()
    BUFFER = 300
    fed = 0
    recieved = 0
    for coreNum in range(cores):
        p = multiprocessing.Process(target=distCullThread, args=(tasks,results,tree,maxCut,resNum))
        jobs.append(p)
    while (fed < BUFFER) and (fed < numberOfPoses):
        poseNum = poseList[fed]
        pep = poseData.getCopiedPeptide(poseNum)
        tasks.put(pep)
        fed+=1
    print("Buffer filled")
    countDisp = Counter.Counter(numberOfPoses,"%d/%d distCull")
    for p in jobs:
        p.start()
    while (recieved < numberOfPoses):
        if (tasks.qsize() < BUFFER):
            if (fed < numberOfPoses):
                poseNum = poseList[fed]
                pep = poseData.getCopiedPeptide(poseNum)
                tasks.put(pep)
                fed+=1
            elif (fed < numberOfPoses + cores):
                tasks.put(None)
                fed+=1
        if not results.empty():
            tup = results.get()
            pNum = tup[0]
            retain = tup[1]
            if (retain):
                poseData.good(pNum)
            else:
                poseData.bad(pNum)
            recieved+=1
            countDisp.disp(recieved)         
    for p in jobs:
        p.join()
    PRINT("Calcs done")
def validityCheck(targetTree,step,poseData):
    PRINT("Beginning Cull Script")
    #minCutoff = step["min"]
    minCutoff = step["max"]
    maxCutoff = step["max"]
    #exceptions = step["exceptions"]
    exceptions = 0
    allAtoms = step.get("allAtoms",True)
    PRINT("min:"+str(minCutoff)+" max:"+str(maxCutoff)+" exceptions allowed for max: "+str(exceptions))
    PRINT("Performing validity on "+str(poseData.getLength())+" poses")
    inRange, noRedundancy = calcDistFromPoses(targetTree,poseData,minCutoff,maxCutoff,exceptions,allAtoms)
    PRINT("\tPoses within distance "+str(inRange))
    return noRedundancy

    
def calcDistFromPoses(tree,data,cMin,cMax,cExcep,allAtoms):
    import RedundancyCheck
    RedundancyCheck.setRecorder(recorder)
    PRINT("Beginning Peptide Culling Calculations")
    CAList = []
    count = 0
    countDisp = Counter.Counter(data.getLength(),"%d/%d distCull with validity")
    for n in range(1,data.getLength()+1):
        countDisp.disp(n)
        pep = data.getPeptide(n)
        #PRINT(str(pep))
        #PRINT("CA "+str(len(pep.CAList)))
        retain = pep.allowedDistance(tree,cMin,cMax,cExcep,allAtoms)
        if (retain):
            count+=1
            data.good(n)
        else:
            data.bad(n)
    PRINT("Culling routine done")
    return count,RedundancyCheck.redundancyCheck(data)
