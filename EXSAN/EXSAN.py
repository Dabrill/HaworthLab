#Stock
import multiprocessing
import os
import shutil
import math
import time
import traceback
import json
import random
import copy
#In-house
import PDBTools
import Log
import Fixvar
import ETable
import HBond
import ZMAT
from ZMAT import zmatName,zmatLoc
import KDTree
import PoseScorer
import Permutations
import BaseWaters
import TMD
import WaterCull
import DistCull
import Counter
def PRINT(string):
    print(string)
    recorder.write(string+"\n")        

def getHelperScriptPath(script):
    return CONSTANTS["universalFolder"]+"/"+CONSTANTS["HelperScripts"][script]["path"]+"/"+CONSTANTS["HelperScripts"][script]["file"]    
def scoreControlPeptide(step):
    seq = step["sequence"]
    startRes = step.get("anchor",None)
    model = PDBTools.readPDBFile(CONSTANTS["commonFolder"]+"/"+CONSTANTS["knownStruc"])[0]

    modelResidues = model.listResidues()
    modelSeq = model.getSequence()
    PRINT("Sequence of full Xray peptide:\n\t"+modelSeq)

    overlaps = model.findOverlaps(seq)
    index = PDBTools.pickOverlap(overlaps,startRes)
    if (index == -1):
        raise Exception("Subsequence "+seq+" not found")
    if (index == -2):
        raise Exception("Multiple overlaps found, and insufficient or illegal information to select one overlap\nPlease specify starting residue, counting up from zero (0)\n\t or specify most N / most C terminal hit")
    if (index == -3):
        raise Exception("Only one overlap found, starting a residue "+str(overlaps[0])+" but did not match request start residue "+str(startRes))
    if (index == -4):
        raise Exception("Unknown error picking start residue from overlaps. Number of overlaps found: "+str(len(overlaps)))
    first = modelResidues[index][0]
    last = modelResidues[index+len(seq)-1][0]            
    PRINT("Sequence found from "+str(first)+" to "+str(last))            
            
                    

    atomCount = 1
    pdbName = seq+".pdb"
    pdbOutF = open(pdbName,'w')
    for a in model.atoms:
        if (a.residueNumber >= first) and (a.residueNumber <= last):
            a.residueNumber = a.residueNumber-first+1
            a.atomNumber = atomCount
            a.chain = "L"
            atomCount+=1
            pdbOutF.write(a.toPDBLine()+"\n")
    pdbOutF.close()
            
    P = Permutations.PeptidePermutations.readPDB(pdbName)
    P.saveFile("All.cpf")
def atomCull(fxvr,plan):
    atom = plan["atoms"]
    if ("strict" in plan) and(plan["strict"]):
        fxvr.atomFreeMultB(atom)
        poseNumRecorder.recordStep(Log.QuantityLog.ATM_CULL,fxvr.numGood(),str(atom))
    else:
        for i in range(0,len(atom)):
            sub = []
            for j in range(i,len(atom)):
                sub.append(atom[j])
            PRINT(str(sub))
            fxvr.atomFreeMultB(sub)
            poseNumRecorder.recordStep(Log.QuantityLog.ATM_CULL,fxvr.numGood(),str(sub))
            PRINT("\t"+str(fxvr.numGood()))
def backboneCull(data,step):
    atoms = step["atoms"]

    dataGood = len(data.goodPoseList(True)[0])
    
    PRINT("\tBefore:"+str(data.fixvar.numGood()))
    data.fixvar.backboneCluster(atoms)
    poseNumRecorder.recordStep(Log.QuantityLog.BB_CULL,data.fixvar.numGood(),str(atoms))
    PRINT("\tAfter:"+str(data.fixvar.numGood()))
    dataGood = len(data.goodPoseList(True)[0])
    data.saveCullResults("BackboneCull.pos")
def atomCullReviseOutputs(poseData):
    poseData.saveCullResults("AtomCull.pos")
def applySelectPosesToFixvar(poseData,file):
    PRINT("Editing fixvar file to include only retained poses")
    good = []
    inf = open(file,'r')
    for l in inf:
        good.append(int(l.replace('\n','')))
    good.sort()
    poseData.whitelist(good)  
    poseData.saveCullResults("WaterCull.pos")

def createPeptidePermutations(step,fxvr,tree):
    print("Preparing threads to calculate pose coordinates")
    zObj = ZMAT.ZMAT(zmatLoc(step))   
    zObj.readReference(CONSTANTS)
    orderMap = zObj.createOrderMap()
    exAtms = zObj.convertToCartesian(fxvr,0,"L")
    exAtms = zObj.arrangeAtoms(exAtms)
    example = PDBTools.Peptide(exAtms,0)
    example.makeDictionary()
    pp = Permutations.PeptidePermutations(example,fxvr.getLength())
    pp.linkToFixvar(fxvr)
    intraCulled = 0
    numberOfPoses = fxvr.getLength()
    countDisp = Counter.Counter(numberOfPoses,"%d/%d 3D-Coord")
    

    cores = min(CONSTANTS["cores"],numberOfPoses)
    tasks = multiprocessing.JoinableQueue()
    results = multiprocessing.Queue()
    consumers = [ ZMAT.ThreeDCoord_Thread(zObj,fxvr,CONSTANTS,tree,tasks,results) for i in range(cores) ]   
    BUFFER = 300
    fed = 0
    recieved = 0
    while (fed < BUFFER) and (fed < numberOfPoses):
        tasks.put((fed+1,fxvr.poses[fed]))
        fed+=1
    for w in consumers:
        w.start()
    print(cores)
    print("Buffer filled. Creating pose coordinates")
    while (recieved < numberOfPoses):
        if (tasks.qsize() < BUFFER):
            if (fed < numberOfPoses):
                tasks.put((fed+1,fxvr.poses[fed]))
                fed+=1
            elif (fed < numberOfPoses + cores):
                tasks.put(None)
                fed+=1
        if not results.empty():
            res = results.get()
            rPoseNum = res[0]
            coordMatrix = res[1]
            poseValid = res[2]
            recieved+=1
            pp.feedPoseIntArray(coordMatrix,rPoseNum)
            if not (poseValid):
                intraCulled+=1
                pp.bad(rPoseNum)
            countDisp.disp(recieved)
    PRINT("Number of intramolecular clashes "+str(intraCulled))
    pp.saveCullResults("IntraCull.pos")
    return pp
def sRemove(file):
    try:
        os.remove(file)
        return True
    except:
        return False
def RunWatGen(currentPoseData):
    import Watgen3
    Watgen3.setRecorder(recorder)
    return Watgen3.RunWatGen(CONSTANTS,currentPoseData)
def FakeWatGen(currentPoseData):
    import Watgen3
    Watgen3.setRecorder(recorder)
    return Watgen3.FakeWatGen(CONSTANTS,currentPoseData)
def distCull(tree,step,data,res):
    cores = min(CONSTANTS["cores"],data.numGood())
    PRINT("Running distance cull with %d cores"%cores)
    DistCull.distCull(tree,step,data,res,cores)
    goodPoses = data.numGood()
    PRINT("\tPoses still valid "+str(goodPoses))
    poseNumRecorder.recordStep(Log.QuantityLog.DIST_CULL,goodPoses)
    data.saveCullResults("DistCull.pos")
def getNextBBGrowStep(plan,n):
    '''
    for i in range(len(plan)):
        print(i)
        print("\t",plan[i]["type"])
        if ("sequence" in plan[i]):
            print("\t",plan[i]["sequence"])
    '''
    baseline = plan[n]
    if ("sequence" not in baseline):
        baseline = priorTMDstep(plan,n)
    baselineSize = len(baseline["sequence"])
    i = n
    while (i < len(plan)):
        if plan[i]["type"] in TMD_Like_Steps:
            if (len(plan[i]["sequence"]) > baselineSize):
                return plan[i]
        i+=1
    return None        
def extendibilityCull(targetTree,step,data):
    if step is None:
        return
    import ExtendabilityCull
    forward = (step["reverse"] != 1)
    ExtendabilityCull.extendibilityCull(CONSTANTS,data,step)
    goodPoses = data.numGood()
    PRINT("\tPoses still valid "+str(goodPoses))
    poseNumRecorder.recordStep(Log.QuantityLog.EXTEND_CULL,goodPoses)
    data.saveCullResults("ExtendCull.pos")
    
def validityCheck(targetTree,step,data):
    noRedundancy = DistCull.validityCheck(targetTree,step,data)
    goodPoses = data.numGood()
    if (noRedundancy):
        poseNumRecorder.recordStep(Log.QuantityLog.DIST_CULL,goodPoses,"OK")
    else:
        poseNumRecorder.recordStep(Log.QuantityLog.DIST_CULL,goodPoses,"Redundancy Found")
    PRINT("\tPoses still valid "+str(goodPoses))    
    data.saveCullResults("DistCull.pos")   
def scorePoses(currentPoseData,hydrogenPoseData,waterData,Zmat,BaseWaters):
    def makeTaskTuple(n):
        pep = currentPoseData.poseBit(n)
        wat = waterData.grabPose(n,True)
        #Pose Number, Pose Peptide, Pose Water
        return (n,pep,wat)
    from ScoreFile import ScoreFile
    PRINT("Running multicore pose scoring routine")
    #fxvr = currentPoseData.fixvar
    
    poseList,includeAll = currentPoseData.goodPoseList()
    #XXXXXXXX
    #poseList = sorted(random.sample(poseList,1000))
    #poseList = poseList[:65000]
    '''
    poseList = []
    for p in range(1,5000):
        poseList.append(p)
    '''
    dataGood = len(currentPoseData.goodPoseList(True)[0])
    #input("CurPoseData %i/%i"%(dataGood,len(poseList)))
    #example = hydrogenPoseData.getPeptide(1)
    example = currentPoseData.getPeptide(1)
    example.makeDictionary()
    numberOfPoses = len(poseList)        
    scoreTable = ETable.EnergyTable.makeBlankTable(currentPoseData.getLength(),"ETable.tsv")
    topFolder = os.path.abspath('')
    alignDump = ScoreFile(poseList,"%s/alignDump.dtf"%topFolder)
    sideDump = ScoreFile(poseList,"%s/sideDump.dtf"%topFolder)
    hbondFile = ScoreFile(poseList,"%s/HBondTable.dtf"%topFolder)
    
    cores = min(CONSTANTS["cores"],numberOfPoses)
    #cores = 1
    tasks = multiprocessing.JoinableQueue()
    results = multiprocessing.Queue()
    hBondInfo = multiprocessing.Queue()
    alignInfo = multiprocessing.Queue()
    sideChainInfo = multiprocessing.Queue()
    consumers = [ PoseScorer.PoseScorer(example,tasks, results, hBondInfo, alignInfo,sideChainInfo, CONSTANTS, Zmat, BaseWaters) for i in range(cores) ]   
    BUFFER = 300
    fed = 0
    recieved = 0
    
    #DEBUG XXX
    #numberOfPoses = 10000
    while (fed < BUFFER) and (fed < numberOfPoses):
        poseNum = poseList[fed]
        tup = makeTaskTuple(poseNum)
        tasks.put(tup)
        fed+=1
    for w in consumers:
        w.start()
    print("Buffer filled. Running scoring")
    start = time.time()
    scoreStart = time.time()
    countDisp = Counter.Counter(numberOfPoses,"%d/%d score")
    while (recieved < numberOfPoses):
        if (tasks.qsize() < BUFFER):
            if (fed < numberOfPoses):
                poseNum = poseList[fed]
                tup = makeTaskTuple(poseNum)
                tasks.put(tup)
                fed+=1
            elif (fed < numberOfPoses + cores):
                tasks.put(None)
                fed+=1
        if not alignInfo.empty():
            alignDump.feedPose(alignInfo.get())
        if not hBondInfo.empty():
            hbondFile.feedPose(hBondInfo.get())
        if not sideChainInfo.empty():
            sideDump.feedPose(sideChainInfo.get())
        if not results.empty():
            res = results.get()
            resPoseNum = res[0]
            values = res[1]
            for cell in values:
                if (cell[1] == "Energy"):
                    currentPoseData.fixvar.setEnergy(resPoseNum,cell[0])
                scoreTable.feedParameter(cell[0],cell[1],resPoseNum)
            recieved+=1
            countDisp.disp(recieved)
            if (recieved % 10000 == 0):
                recorder.write("\t%s     %i\t%2.3f\n"%(Log.timeString(),recieved,time.time()-start))
                start = time.time()
                recorder.copyLog()
    tasks.join()        
    PRINT("Calculations complete")
    PRINT("Scoring took "+str(time.time()-scoreStart)+" seconds total")
    PRINT("num active children:"+str(multiprocessing.active_children()))
    PRINT("Outputing table")
    scoreTable.doneFeeding()
    alignDump.end()
    hbondFile.end()
    sideDump.end()
    PRINT("Done processing")
    currentPoseData.energyTable = scoreTable
    try:
        currentPoseData.fixvar.energyInfo = True
    except:
        PRINT("NO FIXVAR VALUE TO APPLY ENERGY TO")

def makeLMSTableWithoutWatgen(poseData,anchor):
    model = PDBTools.readPDBFile(CONSTANTS["commonFolder"]+"/"+CONSTANTS["knownStruc"])[0]
    model.makeDictionary()
    if (anchor == "C"):
        start = len(model.residues) - 1
    else:
        firstPep = poseData.getPeptide(0)
        overlaps = model.findOverlaps(firstPep)
        start = PDBTools.pickOverlap(overlaps,anchor)

    firstPep = poseData.getPeptide(0)
    #overlaps = model.findOverlaps(firstPep)
    #start = anchor
    missingAtoms = firstPep.compareMissing(model,start)
    if(len(missingAtoms) > 0):
        PRINT("Following atoms could not be located in Xray control:")
        for ma in missingAtoms:
            PRINT("\t"+ma.identifier())

    #print(model.atoms[62].identifier())
    distRef = None
    for a in model.atoms:
        if (a.atomType == "C"):
            distRef = a
    print(distRef.identifier())
    print(firstPep.atoms[0].identifier())
    
    poseList,includeAll = poseData.goodPoseList()
    outf = open("ControlCompare.tsv",'w')
    outf.write("Pose\tLMS\tLMS_BB\tDistance\n")
    for p in poseList:
        if (p % 25000 == 0):
            print(p)
        pep = poseData.getPeptide(p)
        if (start >= 0):
            LMS, LMS_BB = pep.compare(model,start)
            outf.write(str(p)+'\t'+str(LMS)+'\t'+str(LMS_BB)+'\t'+str(distRef.distance(pep.atoms[0]))+'\n')
        else:
            outf.write(str(p)+'\t'+str(start)+'\t'+str(start)+'\n')
    outf.close()


def cutControlTable(cullType):
    cullFilename = cullType+"Cull.pos"
    poseList = Permutations.readGoodPoseFile(cullFilename)
    wholeTable = ETable.EnergyTable.readTableFile("ConTable.tsv")
    PRINT("Control Table loaded")
    subTable = wholeTable.makeSubtable("CulledTable.tsv",poseList)
    subTable.rankTable("Energy","E_Rank")
    subTable.rankTable("RMS","RMS_Rank")
    subTable.rankTable("RMS_BB","RMS_BB_Rank")
    subTable.doneFeeding()
def consensusPose(poseBank):
    def outputConsensus(poseline,tier,coveredFamilies):
        pNum = poseline[0]
        filename = "Consensus_Pose%i_Tier%i.pdb"%(pNum,tier)
        outf = open(filename,'w')
        outf.write(poseBank.getPosePDB(pNum,False))
        outf.close()
        for k in poseline[3]:
            coveredFamilies.add(k)
    SIZE = 100
    CUTOFF = 2.0
    MIN_FAMILY = 19
    ONLY_OUTPUT_CUTOFF = 50
    goodPoses = poseBank.goodPoseList()[0]
    energyTable = []
    SIZE = min(SIZE,len(goodPoses))
    for pNum in goodPoses:
        tup = [pNum,poseBank.fixvar.getEnergy(pNum),None,None]
        energyTable.append(tup)
    energies = sorted(energyTable, key=lambda tup: tup[1])
    bestPoses = energies[:SIZE]
    energies = None
    tablesPreferenceOrder = ["CulledTable.tsv","ConTable.tsv","ETable.tsv"]
    inTable = None
    for filename in tablesPreferenceOrder:
        if (os.path.isfile(filename)):
            PRINT("Adding consensus column to table "+filename)
            inTable = ETable.EnergyTable.readTableFile(filename)
            break
    if inTable is not None:
        inTable.addCols(["Consensus"],[int])
        for i in range(len(bestPoses)):
            entry = bestPoses[i]
            pNum = entry[0]
            pose = poseBank.getCopiedPeptide(pNum)
            entry[2] = pose
            entry[3] = []

    for i in range(len(bestPoses)):
        poseA = bestPoses[i][2]
        for j in range(i+1,len(bestPoses)):
            poseB = bestPoses[j][2]
            #if (poseA is not poseB):
            conScore = poseA.compare(poseB)[1]
            if (conScore <= CUTOFF):
                bestPoses[i][3].append(bestPoses[j][0])
                bestPoses[j][3].append(bestPoses[i][0])
    bestPoses = sorted(bestPoses, key=lambda tup: tup[0])
    '''
    for e in bestPoses:
        print(e[0],e[1],e[2].indexNum,e[3])
    '''

    
    pos = 0
    for i in range(inTable.length()):
        tablePNum = inTable.getParameter("Pose",i)
        if (pos < len(bestPoses)):
            bPNum = bestPoses[pos][0]
        else:
            bPNum = None
        if (bPNum == tablePNum):
            inTable.feedParameter(len(bestPoses[pos][3]),"Consensus",i)
            pos+=1
        else:
            inTable.feedParameter(0,"Consensus",i)

    inTable.writeToFile("ConsensusTable.tsv")
    bestPoses = sorted(bestPoses, key=lambda tup: tup[1])
    bestPoses = sorted(bestPoses, key=lambda tup: len(tup[3]), reverse=True)
    coveredFamilies = set()
    outputConsensus(bestPoses[0],1,coveredFamilies)
    if (len(bestPoses[0][3]) > ONLY_OUTPUT_CUTOFF):
        return
    tier = 2
    for i in range(1,SIZE):
        pose = bestPoses[i]
        if (pose[0] not in coveredFamilies):
            if (len(pose[3]) >= MIN_FAMILY):
                outputConsensus(pose,tier,coveredFamilies)
                tier+=1


def priorGly(p,n):
    i = n-1
    while (i>=0):
        if p[i]["type"] == "TMD":
            if (p[i]["reverse"]):
                if(p[i]["sequence"][0] == "G"):
                    return p[i]
            else:
                seq = p[i]["sequence"]
                if(seq[len(seq)-1] == "G"):
                    return p[i]
        i-=1
    return "__START__"
TMD_Like_Steps = ["TMD","initCTerminal"]
def isBBGrow(p,n):
    length = len(p[n]["sequence"])
    for i in range(n-1,-1,-1):
        if (p[i]["type"] in TMD_Like_Steps):
            nLength = len(p[i]["sequence"])
            if (length > nLength):
                return True
            else:
                return False
    return True
def priorTMDstep(p,n):
    i = n-1
    while (i>=0):
        if p[i]["type"] in TMD_Like_Steps:
            return p[i]
        i-=1
    return "__START__"
def priorSeq(p,n,returnString = True):
    i = n
    while (i>=0):
        if "sequence" in p[i]:
            if returnString:
                return p[i]["sequence"]
            else:
                return p[i]
        i-=1
    return "__START__"
def priorPlanFolder(p,n):
    i = n-1
    while (i>=0):
        if "folder" in p[i]:
            return top+"/"+p[i]["folder"]
        i-=1
    return "__START__"
def priorCullStep(p,n):
    i = n-1
    while (i>=0):
        if p[i]["type"] == "TMD":
            for possible in ("Extend","Intra"):
                cullFilename = possible+"Cull.pos"
                if os.path.isfile(cullFilename):
                    return possible
        if p[i]["type"] == "backboneCull":
            return "Backbone"
        if p[i]["type"] == "waterCull":
            return "Water"
        if p[i]["type"] == "atomCull":
            return "Atom"
        if p[i]["type"] == "distCull" or p[i]["type"] == "validityCheck":
            return "Dist"
        i-=1
    return None
def executePlan(plan):
    global poseNumRecorder
    start = 0
    stop = len(plan)
    PRINT("Preparing prescript variables")
    for n in range(start,len(plan)):
        if "begin" in plan[n]:
            if (plan[n]["begin"]) and (start == 0):
                start = n
        if "end" in plan[n]:
            if plan[n]["end"]:
                if (n > start and n < stop):
                    stop = n
    PRINT("Start point determined")
    commonFolder = CONSTANTS["commonFolder"]
    targetFilename = CONSTANTS["targetProtein"]
    targetPDBFile = "%s/%s"%(commonFolder,targetFilename)
    #radius = CONSTANTS["cullRadiusMeasure"]
    #centerResidue = CONSTANTS["cullRadiusResidue"]
    permLogLoc = CONSTANTS.get("permLogLoc",None)
    #PRINT("Target Alpha Carbon consideration radius: "+str(radius)+" from residue "+str(centerResidue))

    #Creates a searchable tree with all the atoms in the target protein
    targetAll = PDBTools.Peptide(PDBTools.readPDBSelectively(targetPDBFile,CONSTANTS["targetProteinChain"]))#PDBTools.readPDBFile(commonFolder+"/"+targetFilename)[0]
    targetTree = KDTree.KDTree.loadAtomArray(targetAll.atoms)

    PRINT("Ensuring base water file exists")
    if BaseWaters.shouldGenerateBaseWaterFile(CONSTANTS):
        PRINT("Creating base waters")
        baseWaters = BaseWaters.makeBaseWaters(CONSTANTS,targetTree,targetAll)
    baseWaters = BaseWaters.setupBaseWaters(CONSTANTS)
    PRINT("Base Waters Scored")


    fixvarfile = None
    currentPoseData = None
    waterData = None
    if (start > 0):
        PRINT("Prior Quant Record "+str(start))
        seq = priorSeq(plan,start)
        currentFolder = priorPlanFolder(plan,start)
        if (currentFolder != "__START__"):
            poseNumRecorder = Log.QuantityLog(plan[start],seq)
            poseNumRecorder.top = top
            os.chdir(currentFolder)
            PRINT("Loading Fixvar")
            fixvarfile = Fixvar.Fixvar()
            fixvarfile.setBB(ZMAT.ZMAT(zmatLoc(priorTMDstep(plan,start))))
            PRINT("Trying to load CPF")
            try:
                currentPoseData = Permutations.PeptidePermutations.readCPF("all.cpf")
                currentPoseData.linkToFixvar(fixvarfile)
            except:
                PRINT(traceback.format_exc())
                PRINT("Could not load CPF")       
            priorCullType = priorCullStep(plan,start)
            if (priorCullType is not None):
                cullFilename = priorCullType+"Cull.pos"
                print("PCT ",priorCullType)
                PRINT("Loading which poses were culled during last "+priorCullType+" cull")
                currentPoseData.applySavedCull(cullFilename)
            else:
                print("No prior cull")
            #if not "sequence" in plan[start]:
            #poseNumRecorder.setSequence()
            #poseNumRecorder.alignStep(plan[start]):
    else:
        poseNumRecorder = Log.QuantityLog()
        PRINT("New Quant Record")
    for n in range(start,stop):
        if (plan[n]["type"] in TMD_Like_Steps):
            currentFolder = top+"/"+plan[n]["folder"]
            try:
                os.stat(currentFolder)
            except:
                os.mkdir(currentFolder)
            PRINT("\n\n\n")
            poseNumRecorder.setSequence(plan[n]["sequence"])
            recorder.write("********************\n")
            PRINT("Starting Script on Peptide "+plan[n]["sequence"])
            recorder.write("********************\n")
            os.chdir(currentFolder)
            PRINT(Log.timeString())
        if (plan[n]["type"] == "TMD"):
            ZMAT.createZMAT(zmatLoc(plan[n]))
            numNewPoses = TMD.runTMD(CONSTANTS,plan[n],fixvarfile,recorder)
            poseNumRecorder.recordStep(Log.QuantityLog.TMD,numNewPoses)
            PRINT("TMD Calcs finished at "+Log.timeString())
            PRINT("\tPoses created "+str(numNewPoses))
            if isBBGrow(plan,n):
                if (n == 0):
                    fixvarfile = Fixvar.Fixvar()
                PRINT(Log.timeString())
                recorder.copyLog()
                seeker = HBond.PerfectHydrogenBondSeeker(CONSTANTS,plan[n],fixvarfile,targetTree)
                seeker.analyzeAllPoses()
                fixvarfile = seeker.conclude()
                del seeker
                poseNumRecorder.recordStep(Log.QuantityLog.HBSEEK,fixvarfile.numGood())

            else:
                fixvarfile = Fixvar.Fixvar()
                fixvarfile.setBB(ZMAT.ZMAT(zmatLoc(plan[n])))
            
            PRINT("")
            recorder.copyLog()
            PRINT(Log.timeString())
            #makeOnePDB(plan[n],fixvarfile)
            PRINT("Calculating Pose Coordinates")
            #currentPoseData = Permutations.PeptidePermutations.readPDB("makeOnePDB/All.pdb")
            currentPoseData = createPeptidePermutations(plan[n],fixvarfile,targetTree)
            #PRINT("GOOD POSES: "+str(currentPoseData.numGood()))
            #currentPoseData.linkToFixvar(fixvarfile)
            PRINT("Saving dense pose file")

            currentPoseData.saveFile("all.cpf")
            if isBBGrow(plan,n):
                PRINT("Running extendibility cull")
                extendibilityCull(targetTree,getNextBBGrowStep(plan,n),currentPoseData)
  
            shutil.copy(getHelperScriptPath("GUI"),CONSTANTS["HelperScripts"]["GUI"]["file"])
            PRINT("")
            recorder.copyLog()
            PRINT(Log.timeString())
            PRINT("TMD step Complete")
        elif (plan[n]["type"] == "initCTerminal"):
            import Initializer
            Initializer.initCTerminal(plan[n],CONSTANTS,targetAll)
            PRINT("Initial Routine Sucessful")
            fixvarfile = Fixvar.Fixvar()
            PRINT("")
            recorder.copyLog()
            PRINT(Log.timeString())
            PRINT("Calculating Pose Coordinates")
            currentPoseData = createPeptidePermutations(plan[n],fixvarfile,targetTree)
            PRINT("GOOD POSES: "+str(currentPoseData.numGood()))
            #currentPoseData.linkToFixvar(fixvarfile)
            PRINT("Saving dense pose file")

            currentPoseData.saveFile("all.cpf")
  
            shutil.copy(getHelperScriptPath("GUI"),CONSTANTS["HelperScripts"]["GUI"]["file"])
            PRINT("")
            recorder.copyLog()
            PRINT(Log.timeString())
            PRINT("Initiation step Complete")
        else:
            if "folder" in plan[n]:
                currentFolder = top+"/"+plan[n]["folder"]
                os.chdir(currentFolder)
            else:
                currentFolder = priorPlanFolder(plan,n)
                if (currentFolder != "__START__"):
                    os.chdir(currentFolder)
        if (plan[n]["type"] == "seekHB"):
            raise Exception("Classic Hbond seek feature removed")
        if (plan[n]["type"] == "water"):
            recorder.write("********************\n")
            PRINT("Starting water on Peptide "+priorSeq(plan,n))
            recorder.write("********************\n")
            PRINT(Log.timeString())
            recorder.copyLog()
            highRes = PDBTools.lastResidueNum()
            hydrogenPoseData,waterData = RunWatGen(currentPoseData)
            PRINT(Log.timeString())
            os.chdir(currentFolder)
            zzz = ZMAT.ZMAT(zmatLoc(priorSeq(plan,n,False)))
            zzz.readReference(CONSTANTS)
            scorePoses(currentPoseData,hydrogenPoseData,waterData,zzz,baseWaters)
            #RunWatAnalysis(currentPoseData,hydrogenPoseData,waterData,fixvarfile,targetTree)
            PRINT(Log.timeString())
            os.chdir(currentFolder)     
            PRINT("")
            recorder.copyLog()
            os.chdir(currentFolder)
            #applyEnergyToFixvar(currentPoseData)
        if (plan[n]["type"] == "noWater"):
            recorder.write("********************\n")
            PRINT("Faking water on Peptide "+priorSeq(plan,n))
            recorder.write("********************\n")
            PRINT(Log.timeString())
            recorder.copyLog()
            highRes = PDBTools.lastResidueNum()
            hydrogenPoseData,waterData = FakeWatGen(currentPoseData)
            PRINT(Log.timeString())
            os.chdir(currentFolder)
            zzz = ZMAT.ZMAT(zmatLoc(priorSeq(plan,n,False)))
            zzz.readReference(CONSTANTS)
            scorePoses(currentPoseData,hydrogenPoseData,waterData,zzz,baseWaters)
            #RunWatAnalysis(currentPoseData,hydrogenPoseData,waterData,fixvarfile,targetTree)
            PRINT(Log.timeString())
            os.chdir(currentFolder)     
            PRINT("")
            recorder.copyLog()
            os.chdir(currentFolder)
            #applyEnergyToFixvar(currentPoseData)
        if (plan[n]["type"] == "analysis"):
            recorder.write("********************\n")
            PRINT("Starting only analysis on Peptide "+priorSeq(plan,n))
            recorder.write("********************\n")
            if (waterData is None):
                PRINT("NEED TO LOAD CURRENT POSE DATA")
                #currentPoseData = Permutations.PeptidePermutations.readCPF("All.cpf")
                #currentPoseData.linkToFixvar(fixvarfile)
                hydrogenPoseData  = Permutations.PeptidePermutations.readCPF("AllH.cpf")
                waterData = Permutations.WaterPermutations.readCWF("Water.cwf")
                PRINT("FILES LOADED")
            PRINT(Log.timeString())
            os.chdir(currentFolder)           
            #RunWatAnalysis(currentPoseData,hydrogenPoseData,waterData,fixvarfile,targetTree)
            zzz = ZMAT.ZMAT(zmatLoc(priorSeq(plan,n,False)))
            zzz.readReference(CONSTANTS)
            scorePoses(currentPoseData,hydrogenPoseData,waterData,zzz,baseWaters)
            PRINT(Log.timeString())
            os.chdir(currentFolder)     
            PRINT("")
            recorder.copyLog()
            os.chdir(currentFolder)
            #applyEnergyToFixvar(currentPoseData)
        if (plan[n]["type"] == "atomCull"):
            if not fixvarfile.energyInfo:
                if (os.path.isfile("ETable.tsv")):
                    currentPoseData.applyEnergyToFixvar()
            recorder.write("********************\n")
            PRINT("Culling poses by side-chain redundancy")
            recorder.write("********************\n")
            atomCull(fixvarfile,plan[n])
            if (n < len(plan) - 1):
                if not (plan[n+1]["type"] == "atomCull"):
                    PRINT("Saving output")
                    atomCullReviseOutputs(currentPoseData)
            PRINT("Done with cull")
            PRINT("\tRetained poses: "+str(fixvarfile.numGood()))
        if (plan[n]["type"] == "backboneCull"):
            if not fixvarfile.energyInfo:
                if (os.path.isfile("ETable.tsv")):
                    currentPoseData.applyEnergyToFixvar()
                else:
                    raise Exception("Must have energy data to run this cull")
            recorder.write("********************\n")
            PRINT("Culling poses by backbone similiarity")
            recorder.write("********************\n")
            backboneCull(currentPoseData,plan[n])
            PRINT("Done with cull")
            PRINT("\tRetained poses: "+str(fixvarfile.numGood()))           
        if (plan[n]["type"] == "ControlCompare"):
            recorder.write("********************\n")
            PRINT("Producing ordered list of resemblance to known structure\nNot using energy data")
            recorder.write("********************\n")
            anchor = plan[n].get("anchor",None)
            makeLMSTableWithoutWatgen(currentPoseData,anchor)
        if (plan[n]["type"] == "KnownCompare"):
            recorder.write("********************\n")
            PRINT("Producing ordered list of resemblance to known structure")
            recorder.write("********************\n")
            anchor = plan[n].get("anchor",None)
            ETable.extendETable(CONSTANTS,currentPoseData,anchor)
            priorCullType = priorCullStep(plan,n)
            if priorCullType is not None:
                PRINT("PriorCullType "+str(priorCullType))
                PRINT("Producing table of poses that still remain")
                cutControlTable(priorCullType)
        if (plan[n]["type"] == "waterCull"):
            if not fixvarfile.energyInfo:
                if (os.path.isfile("ETable.tsv")):
                    currentPoseData.applyEnergyToFixvar()
                else:
                    raise Exception("Cannot perform water cull without WATGEN and scoring")
            recorder.write("********************\n")
            PRINT("Culling peptide by water energies")
            recorder.write("********************\n")            
            if ((plan[n]["method"] == "50P_Energy") or (plan[n]["method"] == "fracMaxEnergy")):
                frac = plan[n]["threshold"]
                PRINT("Threshold of max method with threshold "+str(frac))
                WaterCull.thresholdOfMaxEnergy(plan[n],currentPoseData)
                PRINT("\tNumber of poses retained:"+str(currentPoseData.numGood()))           
                poseNumRecorder.recordStep(Log.QuantityLog.WAT_CULL,currentPoseData.numGood(),"Within "+str(frac*100)+"%")
            if (plan[n]["method"] == "top"):
                topNum = plan[n]["liminalNumber"]
                WaterCull.topWaters(plan[n],currentPoseData)
                PRINT("Top %i poses in terms of water energy"%topNum)
                poseNumRecorder.recordStep(Log.QuantityLog.WAT_CULL,currentPoseData.numGood(),"Top %i"%topNum)
            if (plan[n]["method"] == "rubber"):
                WaterCull.rubberTopWaters(plan[n],currentPoseData)
                PRINT("Performing water cull with \"rubber\" retension")
                poseNumRecorder.recordStep(Log.QuantityLog.WAT_CULL,currentPoseData.numGood(),"Top "+str(currentPoseData.numGood()))
            if (plan[n]["method"] == "stochastic"):
                WaterCull.stochasticWaterCull(plan[n],currentPoseData)
                PRINT("Performing stochastic cull")
                poseNumRecorder.recordStep(Log.QuantityLog.WAT_CULL,currentPoseData.numGood(),"Stochastic")
            if (plan[n]["method"] == "delta"):
                PRINT("Delta Cull")
                delta = WaterCull.getDeltaFrac(top,plan,n)
                PRINT("Culling poses at %3.1f%% Threshold"%(delta*100))
                fakeStep = {"threshold":delta,"relative":True}
                WaterCull.thresholdOfMaxEnergy(fakeStep,currentPoseData)
                PRINT("\tNumber of poses retained:"+str(currentPoseData.numGood()))
                poseNumRecorder.recordStep(Log.QuantityLog.WAT_CULL,currentPoseData.numGood(),"Within %3.1f%%"%(delta*100))
        if (plan[n]["type"] == "distCull"):
            res = plan[n].get("residue",None)
            if res is None:
                lastStep = priorGly(plan,n)
                atmNum = lastStep["variableTorsions"][0]["atom"]
                res = ZMAT.residueNumber(zmatLoc(lastStep),atmNum)
            PRINT("Last res: "+str(res))
            distCull(targetTree,plan[n],currentPoseData,res)
        if (plan[n]["type"] == "overlayCull"):
            if not fixvarfile.energyInfo:
                if (os.path.isfile("ETable.tsv")):
                    currentPoseData.applyEnergyToFixvar()
                else:
                    raise Exception("Cannot perform water cull without WATGEN and scoring")
            recorder.write("********************\n")
            PRINT("Overlay Cull")
            recorder.write("********************\n")
            import OverlayCull
            OverlayCull.findOverlays(currentPoseData)
        if (plan[n]["type"] == "validityCheck"):
            recorder.write("********************\n")
            PRINT("Culling peptide by distance from protein")
            recorder.write("********************\n")
            PRINT(Log.timeString())
            recorder.copyLog()
            validityCheck(targetTree,plan[n],currentPoseData)
            PRINT(Log.timeString())
            recorder.copyLog()
        if (plan[n]["type"] == "gate"):
            recorder.write("********************\n")
            PRINT("Applying logic gate")
            if (plan[n]["method"] == "fewerPoses"):
                PRINT("FEATURE NOT YET IMPLEMENTED")
            recorder.write("********************\n")
        if (plan[n]["type"] == "ScoreControl"):
            recorder.write("********************\n")
            PRINT("Scoring segment of XRay Control")
            scSeq = plan[n]["sequence"]
            scSeqString = "\tSequence: "+scSeq 
            scStartRes = plan[n].get("anchor",None)
            if scStartRes is not None:
                if isinstance(scStartRes,int): 
                    scSeqString+="\t starting at residue "+str(scStartRes)
                elif isinstance(scStartRes,str):
                    if (scStartRes.upper() == "N"):
                        scSeqString+="\t starting with the most N-terminal match"
                    if (scStartRes.upper() == "C"):
                        scSeqString+="\t starting with the most C-terminal match"
            PRINT(scSeqString)
            recorder.write("********************\n")
            currentFolder = top+"/"+scSeq
            try:
                os.stat(currentFolder)
            except:
                os.mkdir(currentFolder)
            os.chdir(currentFolder)
            scoreControlPeptide(plan[n])
            shutil.copy(getHelperScriptPath("GUI"),CONSTANTS["HelperScripts"]["GUI"]["file"])
            shutil.copy(commonFolder+"/out.pdb",currentFolder+"/out.pdb")
            currentPoseData = Permutations.PeptidePermutations.readCPF("all.cpf")

            hydrogenPoseData,waterData = RunWatGen(currentPoseData)
            zzz = ZMAT.ZMAT(zmatLoc(priorSeq(plan,n,False)))
            zzz.readReference(CONSTANTS)
            currentPoseData.fixvar = Fixvar.Fixvar(file = None)
            currentPoseData.fixvar.poses.append(Fixvar.FixvarPoses([],1))
            scorePoses(currentPoseData,hydrogenPoseData,waterData,zzz,baseWaters)
        if (plan[n]["type"] == "CollectControl"):
            os.chdir(top)
            headerWriten = False
            xrayGrandTable = open("XrayGrandTable.tsv",'w')
            for xrayStep in plan:
                if (xrayStep["type"] == "ScoreControl"):
                    xrayInf = open("%s/%s/ETable.tsv"%(top,xrayStep["sequence"]),'r')
                    header = xrayInf.readline()
                    if not headerWriten:
                        xrayGrandTable.write(header)
                        headerWriten = True
                    xrayData = xrayInf.readline().split("\t")
                    xrayData[0] = xrayStep["sequence"]
                    xrayLine = "\t".join(xrayData)
                    xrayGrandTable.write(xrayLine)
            xrayGrandTable.close()
                
        if (plan[n]["type"] == "Consensus"):
            if not fixvarfile.energyInfo:
                if (os.path.isfile("ETable.tsv")):
                    currentPoseData.applyEnergyToFixvar()
            consensusPose(currentPoseData)
    PRINT("Script Complete")
    sRemove(top+"/log_current.txt")
    timeDifference = int(time.time()-Log.startTime)
    PRINT("Time to completetion: "+Log.secondsToHMS(timeDifference))
    PRINT("Time to completetion (seconds): "+str(timeDifference))
    PRINT("\n\nPose Counts:\n"+poseNumRecorder.dump())
    recorder.close()
    poseNumRecorder.close() 

def Main():
    print("Booting output")
    logFilename = "log_"+time.strftime("%m.%d.%y_%H.%M")+".txt"
    global recorder
    recorder = Log.Log(logFilename)
    WaterCull.setRecorder(recorder)
    ETable.setRecorder(recorder)
    BaseWaters.setRecorder(recorder)
    HBond.setRecorder(recorder)
    ZMAT.setRecorder(recorder)
    PRINT("Welcome to EXplicit Solvent ANchored docking program")
    poseNumRecorder = None
    global top
    top = os.path.abspath('')
    PRINT("Loading Constants")
    global CONSTANTS
    with open('constants.json') as data_file:
        CONSTANTS = json.load(data_file)
    permLogLoc = CONSTANTS.get("permLogLoc",None)
    ZMAT.setCommonFolder(CONSTANTS)
    '''
    tpChains = []
    for l in CONSTANTS["targetProteinChain"]:
        tpChains.append(l)
    '''
    PRINT("EXSAN will now look into plans")
    with open('plan.json') as data_file:    
        planFile = json.load(data_file)
    plan = planFile["plan"]
    from MakeZMAT import CTerminalChar
    for p in plan:
        if (p["type"] in TMD_Like_Steps+["ScoreControl"]):
            p["zsequence"] = p["sequence"]
            p["sequence"] = p["sequence"].replace(CTerminalChar,"")
            ZMAT.zmatObj(p)
            
    PRINT("Plan File Read")
    try:
        executePlan(plan)
    except Exception:
        if poseNumRecorder is not None and poseNumRecorder.open:
            poseNumRecorder.close()
        PRINT("")
        PRINT("")
        PRINT("")
        PRINT("Serious Error Occured.")
        try:
            currentFolder = os.path.abspath('')
            PRINT(currentFolder)
            if ("MAKEPDB" in currentFolder.upper()):
                ind = currentFolder.upper().index("MAKEPDB")
                currentFolder = currentFolder[:ind]            
            TMD.removeWithKey(currentFolder,".zmat")
        except:
            y = 4
        errorMessage = traceback.format_exc()
        PRINT(errorMessage)
        recorder.close()
        if permLogLoc is not None:
            shutil.copy(top+"/"+logFilename,permLogLoc+"/"+logFilename)
        time.sleep(120)
        return errorMessage
    if permLogLoc is not None:
        shutil.copy(top+"/"+logFilename,permLogLoc+"/"+logFilename)
    return ""
if __name__ == "__main__":
    Main()
