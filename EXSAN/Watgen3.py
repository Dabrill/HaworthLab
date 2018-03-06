import threading
import subprocess
import os
import sys
import PDBTools
import Permutations
import Log
import Counter
import time
recorder = None
def PRINT(string):
    global recorder
    print(string)
    if recorder is not None:
        recorder.write(string+"\n")
def setRecorder(obj):
    global recorder
    recorder = obj

class WatThread(threading.Thread):
    numThread = 0
    count = 0
    pingStart = None
    progressPingHappening = False
    #errorFile = open("WGTerror.txt",'w')
    def __init__(me,CONSTANTS,assignment,topFolder,thread,target,inputMatrix,poseObject,waterObject):
        threading.Thread.__init__(me)
        me.CONSTANTS = CONSTANTS
        me.assignment = assignment
        me.topFolder = topFolder
        WatThread.numThread+=1
        me.thread = thread
        me.out = poseObject
        me.outW = waterObject
        me.currentPoseData = inputMatrix
        me.target = target
        if  not (me.thread == WatThread.numThread):
            PRINT("Mismatch in number of threads")
    def closeRecord(me):
        WatThread.errorFile.close()
    def run(me):
        wgBatchJar = me.CONSTANTS["universalFolder"]+"/"+me.CONSTANTS["WGBatch"]
        thstr = str(me.thread)
        BATCH_SIZE = 100
        tpChain = me.CONSTANTS["targetProteinChain"]
        for i in range(0,len(me.assignment),BATCH_SIZE):
            try:
                ball = ""
                for j in range(len(me.target)):
                    ball+=me.target[j]
                ball+="__T__"
                ball+=me.currentPoseData.template.fromTemplate().makePDBString(False)
                ball+="__L__"
                #print(ball)
                ball = ball.encode()
                for j in range(BATCH_SIZE):
                    if (i + j < len(me.assignment)):
                        number = me.assignment[i+j]
                        ball+=me.currentPoseData.poseBit(number)
                p = subprocess.Popen(["java","-jar",wgBatchJar,"L","7","0",str(BATCH_SIZE)], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                back = "Failure"
                backError = ""
                backTuple = ("Failure","Failure")
                try:
                    backTuple = p.communicate(ball)
                    back = backTuple[0].decode('utf-8')
                    backError = backTuple[1].decode('utf-8')
                except:
                    PRINT(traceback.format_exc())
                    PRINT("Bad things [%d:%d]"%(me.thread,number))
                    PRINT("Number of threads active "+str(threading.activeCount()))
                    PRINT("<"+back+">")
                    PRINT("["+backError+"]")
                    PRINT("\""+str(backTuple)+"\"")
                    #while(True):
                        #pass
                finally:
                    p.wait()
                if (len(backError) > 0):
                    PRINT("[JAVA]:"+backError)
                back = back.replace("\r","")
                #PRINT("Pre-eval")
                back = back.split("__POSE__")
                '''
                for j in range(len(back)):
                    PRINT(str(j)+"-"+str(me.assignment[i+j]))
                    cutted = back[j].split("\n")
                    PRINT(str(len(cutted)))
                    for k in range(0,len(cutted)):
                        PRINT(cutted[k])
                '''  
                for j in range(BATCH_SIZE):
                    if (i + j < len(me.assignment)):
                        #PRINT(str(j)+"-"+str(me.assignment[i+j]))
                        number = me.assignment[i+j]
                        poseResults = retainWatgenOutput(tpChain,back[j],number)
                        me.out.feedPose(poseResults)
                        me.outW.feedPose(poseResults)
                        #PING Process
                        WatThread.count=WatThread.count+1
                        if (WatThread.count % 10000 == 0):
                            if not WatThread.progressPingHappening:
                                WatThread.progressPingHappening = True
                                recorder.write("\t%s     %i\t%2.3f\n"%(Log.timeString(),WatThread.count,time.time()-WatThread.pingStart))
                                WatThread.pingStart = time.time()
                                #PRINT("Number of threads active "+str(threading.activeCount()))
                                if (recorder is not None):
                                    recorder.copyLog()
                                WatThread.progressPingHappening = False
                        
                WatThread.counter.disp(WatThread.count)
            except Exception:
                PRINT("["+str(me.thread)+":"+str(number)+"]")
                PRINT(traceback.format_exc())
        #os.remove(inName)
def FakeWatGen(CONSTANTS,currentPoseData):
    def addHSlots(pep):
        for r in pep.residues:
            if (r["C"].residueType == "LYS"):
                for i in range(3):
                    HZ = r["C"].copy()
                    key = "HZ%i"%(i+1)
                    HZ.atomType = key
                    r[key] = HZ
                    pep.atoms.append(HZ)
            if (r["C"].residueType == "ARG"):
                HNE = r["C"].copy()
                HNE.atomType = "HNE"
                r["HNE"] = HNE
                pep.atoms.append(HNE)
                for a in range(1,3):
                    for b in range(1,3):
                        HH = r["C"].copy()
                        key = "%iHH%i"%(b,a)
                        HH.atomType = key
                        r[key] = HH
                        pep.atoms.append(HH)
    def relocateH(pep):
        for r in pep.residues:
            if (r["C"].residueType == "LYS"):
                NZ = r["NZ"]
                CE = r["CE"]
                CD = r["CD"]
                for i in range(3):
                    key = "HZ%i"%(i+1)
                    HZ = r[key]
                    HZ.setAtomLocationByZMAT(NZ,CE,CD,1.0,120,60 + 120*i)
            if (r["C"].residueType == "ARG"):
                NE = r["NE"]
                HNE = r["HNE"]
                CD = r["CD"]
                CZ = r["CZ"]
                HNE.setAtomLocationByZMAT(NE,CD,CZ,1.0,120,-180)
                for a in range(1,3):
                    stemKey = "NH%i"%a
                    stem = r[stemKey]
                    for b in range(2):
                        key = "%iHH%i"%(b+1,a)
                        H = r[key]
                        H.setAtomLocationByZMAT(stem,CZ,NE,1.0,120,0+b*180)
                
        
    topFolder = os.path.abspath('')
    examplePose = currentPoseData.getPeptide(1).copy()
    examplePose.makeDictionary()
    addHSlots(examplePose)
    examplePose.atoms = sorted(examplePose.atoms,key = lambda x: x.residueNumber)
    hydrogenPoseData = Permutations.PeptidePermutations(examplePose,currentPoseData.getLength())
    waterData = Permutations.WaterPermutations(currentPoseData.getLength())
    poseList,includeAll = currentPoseData.goodPoseList()
    cnt = Counter.Counter(currentPoseData.getLength(),"%d/%d fakeWatgen")
    for poseNum in poseList:
        newPeptide = currentPoseData.getPeptide(poseNum)
        for natm in newPeptide.atoms:
            atm = examplePose.getSpecificAtom(natm.residueNumber,natm.atomType)
            atm.location = natm.location
        relocateH(examplePose)
        
        poseBytes = bytearray()
        for i in range(len(examplePose.residues)):
            res = examplePose.residues[i]
            for atmType in res:
                #print(i+examplePose.dictOffset,atmType)
                atm = res[atmType]
                ar = Permutations.PeptidePermutations.pdbToByteCoords(atm.toPDBLine())
                poseBytes+=ar
        wgData = (poseBytes,bytearray(),poseNum)
        waterData.feedPose(wgData)
        hydrogenPoseData.feedPose(wgData)
        cnt.disp(poseNum)
    PRINT("Saving Dense Files")
    PRINT("\tPoses")
    hydrogenPoseData.saveFile(topFolder+"/AllH.cpf")
    PRINT("\tWater")
    waterData.saveFile(topFolder+"/Water.cwf")
    PRINT("\tDone")
    return hydrogenPoseData,waterData        
    
def RunWatGen(CONSTANTS,currentPoseData):
    #fxvr = currentPoseData.fixvar
    PRINT("Starting WatGen Script")
    topFolder = os.path.abspath('')
    allFolder = topFolder+"/makeOnePdb"
    poseList,includeAll = currentPoseData.goodPoseList()
    
    outPDB = PDBTools.fileToStringArray("out.pdb")
    TargetProteinWithHydrogen,FirstPose = probeWatgen(CONSTANTS,outPDB,currentPoseData)
    ##############################
    #poseList = list(range(1,501))
            
    posesNum = len(poseList)
    PRINT("Number of poses for WatGen: "+str(posesNum))
    #baseLoc = [int(FirstPose.atoms[0].location.dims[0] * 1000),int(FirstPose.atoms[0].location.dims[1] * 1000),int(FirstPose.atoms[0].location.dims[2] * 1000)]
    hydrogenPoseData = Permutations.PeptidePermutations(FirstPose,currentPoseData.getLength())
    waterData = Permutations.WaterPermutations(currentPoseData.getLength())   
    
    threadCount = CONSTANTS["cores"]
    WatThread.count = 0
    WatThread.pingStart = time.time()
    WatThread.counter = Counter.Counter(posesNum,"%d/%d watgen")
    PRINT("Running Watgen")
    WatThread.numThread = 0        
    threadPoseAssignment = []
    for a in range(0,threadCount):
        threadPoseAssignment.append([])
    for i in range(len(poseList)):
        a = i % threadCount
        threadPoseAssignment[a].append(poseList[i])
    #Start Loop
    threadLock = threading.Lock()
    threads = []
    for i in range(0,threadCount):
        threads.append(WatThread(CONSTANTS,threadPoseAssignment[i], topFolder, i+1,outPDB,currentPoseData,hydrogenPoseData,waterData))
        threads[i].start()
    for t in threads:
        t.join()
    sys.stdout.flush()
    sys.stdin.flush()
    PRINT("Watgen Complete "+Log.timeString())
    if (recorder is not None):
        recorder.copyLog()
    PRINT("Saving Dense Files")
    PRINT("\tPoses")
    hydrogenPoseData.saveFile(topFolder+"/AllH.cpf")
    PRINT("\tWater")
    waterData.saveFile(topFolder+"/Water.cwf")
    PRINT("\tDone")
    return hydrogenPoseData,waterData

def probeWatgen(CONSTANTS,outPDB,currentPoseData):
    PRINT("Getting Hydrogen-filled TP Structure")
    p = subprocess.Popen(["java","-jar",CONSTANTS["universalFolder"]+"/"+CONSTANTS["WGProbe"],"L","7","0"], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    #PRINT(str(currentPoseData.firstIndex()))
    #PRINT("First "+str(currentPoseData.firstIndex()))
    #PRINT(currentPoseData.getPoseStringArray(currentPoseData.firstIndex()))
    #sentPose = currentPoseData.getPeptide(currentPoseData.firstIndex())
    ball = combineArraysIntoString(outPDB,currentPoseData.getPoseStringArray(currentPoseData.firstIndex()))+"__DONE__\r\n"
    #PRINT(ball)
    back = p.communicate(ball.encode())[0].decode('utf-8')
    back = back.replace("\r","")
    #PRINT("_______________________________")
    #PRINT(">"+str(back))
    TargetProteinWithHydrogen = retainWatgenElement(back,"A")
    FirstPose = PDBTools.Peptide.fromStringArray(retainWatgenElement(back,"L"),1)


    '''
    oxt = sentPose.getTerminalOxygen()
    if oxt is not None:
        FirstPose.atoms.append(oxt)
    '''
    return TargetProteinWithHydrogen,FirstPose
def combineArraysIntoString(in1, in2):
    ret = ""
    for line in in1:
        ret+=line
    #out.write("TER\n")
    for line in in2:
        ret+=line
    return ret
def retainWatgenElement(string,element):
    entries = string.split("\n")
    ret = []
    for line in entries:
        if (element == line[21:22]):
            ret.append(line+"\n")
    return ret
def retainWatgenOutput(targetChain,string,index):
    entries = string.split("\n")
    water = bytearray()
    ligand = bytearray()
    for i in range(len(entries)):
        line = entries[i]
        if (not (targetChain == line[21:22]))  and (len(line) > 3):
            #ret+=line+"\n"
            try:
                ar = Permutations.PeptidePermutations.pdbToByteCoords(line)
            except:
                #PRINT(line[13:15]+"\n")
                if (line[13:14] == "H"):
                    if (line[13:15] == "H1"):
                        t = entries[i-1]
                    if (line[13:15] == "H2"):
                        t = entries[i-2]
                    ar = Permutations.PeptidePermutations.pdbToByteCoords(t)
                else:
                    PRINT("PARSE ERROR")
                    PRINT(line)
                    PRINT("\t"+line[31:38])
                    PRINT("\t"+line[38:46])
                    PRINT("\t"+line[46:54])
                    raise Exception("Parse Error")
                
            if (line[21:22] == "L"):
                ligand = ligand + ar
            if (line[21:22] == "I") or (line[21:22] == "W"):
                water = water + ar
    return ligand,water,index
