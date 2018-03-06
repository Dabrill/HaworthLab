import PDBTools
import BaseWaters
import Permutations
import multiprocessing
import Counter
from KDTree import KDTree
OData = (1.22,119.99,-180.00)

rNData = (1.46,109.42,-179.98)
rCAData = (1.51,120.05,-179.98)
rCData = (1.35,120.03,-180.00)


fNData = (1.35,120.05,-180.00)
fCAData = (1.46,120.03,-179.98)
fCData = (1.51,109.40,-180.00)
class ExtendibilityCuller(multiprocessing.Process):
    def __init__(me,cons,forward,examplePose,phi,psi,task_queue,result_queue):
        multiprocessing.Process.__init__(me)
        try:
            me.target = PDBTools.readPDBFile("out.pdb")[0]
        except:
            me.target = PDBTools.readPDBFile("%s/out.pdb"%cons["commonFolder"])[0]
        me.targetTree = KDTree.loadAtomArray(me.target.atoms)
        me.clashLimit = cons["clashLimit"]
        me.forward = forward
        me.task_queue = task_queue
        me.result_queue = result_queue

        me.phi = phi
        me.psi = psi

        me.examplePose = examplePose
        me.iN = None
        me.iCA = None
        me.iC = None
        highestRes = max(examplePose.atoms,key = lambda x: x.residueNumber).residueNumber
        for i in range(len(examplePose.atoms)):
            atm =  examplePose.atoms[i]
            if (forward and (atm.residueNumber == highestRes)) or ((not forward) and (atm.residueNumber == 1)):
                if (atm.atomType == "N"):
                    me.iN = i
                if (atm.atomType == "C"):
                    me.iC = i
                if (atm.atomType == "CA"):
                    me.iCA = i
                if (atm.atomType == "O"):
                    me.iO = i
    def run(self):
        while(True):
            next_task = self.task_queue.get()
            if next_task is None:
                # Poison pill means shutdown
                #print ('%s: Exiting' % self.name)
                self.task_queue.task_done()
                break
            poseNum = next_task[0]
            poseBytes = next_task[1]
            pose = Permutations.PeptidePermutations.applyBytesToPeptide(self.examplePose,poseBytes)
            result = self.isExtendable(pose)
            self.result_queue.put((poseNum,result))
        return   
    def isExtendable(me,ligand):
        N = ligand.atoms[me.iN]
        CA = ligand.atoms[me.iCA]
        C = ligand.atoms[me.iC]
        O = ligand.atoms[me.iO]
        foundGood = False
        
        if (me.forward):
            for A in range(me.psi["start"],me.psi["stop"]+1,me.psi["step"]):
                N0 = C.location.zmatToCartesian(C.location,CA.location,N.location,fNData[0],fNData[1],A)
                O.location = C.location.zmatToCartesian(C.location,CA.location,N0,OData[0],OData[1],OData[2])
                if me.hasClash([N0,O.location]):
                    #print("N Clash")
                    continue
                for B in range(me.phi["start"],me.phi["stop"]+1,me.phi["step"]):
                    CA0 = C.location.zmatToCartesian(N0,C.location,CA.location,fCAData[0],fCAData[1],fCAData[2])
                    C0 = C.location.zmatToCartesian(CA0,N0,C.location,fCData[0],fCData[1],B)
                    if me.hasClash([CA0,C0]):
                        #print("C Clash")
                        continue
                    return True                                      
        else:
            for A in range(me.phi["start"],me.phi["stop"]+1,me.phi["step"]):
               C0 = C.location.zmatToCartesian(N.location,CA.location,C.location,rCData[0],rCData[1],A)
               CA0 = C.location.zmatToCartesian(C0,N.location,CA.location,rCAData[0],rCAData[1],rCAData[2])
               if me.hasClash([C0,CA0]):
                   continue
               for B in range(me.psi["start"],me.psi["stop"]+1,me.psi["step"]):
                   #print(A,B)
                   N0 = C.location.zmatToCartesian(CA0,C0,N.location,rNData[0],rNData[1],B)
                   O0 = C.location.zmatToCartesian(C0,CA0,N.location,OData[0],OData[1],OData[2])
                   if me.hasClash([N0,O0]):
                       continue
                   return True
               '''
               print(A,B)
               for atm in [C0,CA0,N0,O0]:
                   neigh, dist = targetTree.nearestNeighbor(atm.location,True)
                   print(atm,neigh,dist)
               '''
        return False
    def hasClash(me,atoms):
       for atm in atoms:
           neigh, dist = me.targetTree.nearestNeighbor(atm,True)
           if (dist < me.clashLimit):
               return True
       return False
def extendibilityCull(cons,curPoses,step):
    import multiprocessing
    import ZMAT
    forward = (step["reverse"] != 1)
    #ZMAT.setCommonFolder(cons)
    zmat = ZMAT.zmatObj(step)
    zmat.readReference(cons)
    for i in range(len(step["variableTorsions"])):
        tors = step["variableTorsions"][i]
        rotAtmNum = tors["atom"]
        rotAtm = zmat.zAtoms[rotAtmNum-4]            
        if (tors["type"] == 0):  
            if (rotAtm.atomType == "C"):
                phi = tors
            if (rotAtm.atomType == "N"):
                psi = tors
    #print(phi,psi)
    example = curPoses.getPeptide(1)

    #input("Forward %s"%forward)

    poseList,includeAll = curPoses.goodPoseList()
    numberOfPoses = len(poseList)
    cores = min(cons["cores"],numberOfPoses)
    tasks = multiprocessing.JoinableQueue()
    results = multiprocessing.Queue()

    #input("A")
    consumers = [ ExtendibilityCuller(cons,forward,example.copy(),phi,psi,tasks,results) for i in range(cores) ]
    BUFFER = 300
    fed = 0
    recieved = 0
    good = 0
    bad = 0
    while (fed < BUFFER) and (fed < numberOfPoses):
        poseNum = poseList[fed]
        pose = curPoses.poseBit(poseNum)
        tasks.put((poseNum,pose))
        fed+=1
    for w in consumers:
        w.start()
    countDisp = Counter.Counter(numberOfPoses,"%d/%d extendibility")
    while (recieved < numberOfPoses):
        if (tasks.qsize() < BUFFER):
            if (fed < numberOfPoses):
                poseNum = poseList[fed]
                pose = curPoses.poseBit(poseNum)
                tasks.put((poseNum,pose))
                fed+=1            
            elif (fed < numberOfPoses + cores):
                tasks.put(None)
                fed+=1
        if not results.empty():
            res = results.get()
            resPoseNum = res[0]
            isGoodPose = res[1]
            if isGoodPose:
                curPoses.good(resPoseNum)
                good+=1
            else:
                curPoses.bad(resPoseNum)
                bad+=1
            recieved+=1
            countDisp.disp(recieved)
    print("Done with %i good and %i bad "%(good,bad))
'''   
if __name__ == "__main__":
    import time
    #pp = Permutations.PeptidePermutations.readCPF("%s/all.cpf"%Folder)
    pp = Permutations.PeptidePermutations.readCPF("all.cpf")
    fakeDict = {"clashLimit":1.7,"cores":10}
    start = time.time()
    A = extendibilityCullSingle(fakeDict,pp)
    print("Single ",time.time()-start)
    start = time.time()
    B = extendibilityCull(fakeDict,pp)
    print("Multi ",time.time()-start)
    input("Multi Complete")
    input("Routine complete")
    for key in range(1,1001):
        try:
            a = A[key]
            b = B[key]
            print(key,b)
            if (a!=b):
                print("Mismatch at key %i"%key)
        except:
            input("Error with key %i"%key)
    input("DOne")
if False:# __name__ == "__mainByte__":
    curPoses = Permutations.PeptidePermutations.readCPF("all.cpf")
    five = curPoses.getCopiedPeptide(5)
    eight = curPoses.getCopiedPeptide(800)
    eightBytes = curPoses.poseBit(800)
    #print(eightBytes)
    print(five)
    input("")
    print(eight)
    input("")
    curPoses.applyBytesToPeptide(five,eightBytes)
    print(five)
    input("")    
'''


