import PDBTools
import ZMAT
SidechainData = None

CTerminalChar = "#"
def loadJSON():
    global SidechainData
    import os
    import json
    scriptLocation = os.path.realpath(__file__)
    scriptPath = os.path.dirname(scriptLocation)
    jsonLocation = "%s/%s"%(scriptPath,'ZmatSidechains.json')
    with open(jsonLocation) as data_file:    
        SidechainData = json.load(data_file)
loadJSON()

def SeqToZMAT(taggedSequence,centerResidue,outputFilename = None):
    def popAtom(resNum,atmName):
        z = ZMAT.ZMAT_Atom.blank()
        z.resNum = resNum
        z.resType = seqResidues[resNum]
        z.atomType = atmName
        return z
    def applyConnect(zAtm,data,resNum):
        con = data["Connect"]
        zAtm.distAtom = zAtomDict[(resNum,con[0])].fixvarNum
        zAtm.angAtom = zAtomDict[(resNum,con[1])].fixvarNum
        zAtm.torAtom = zAtomDict[(resNum,con[2])].fixvarNum
        applyGeo(zAtm,data)
    def applyGeo(zAtm,data):
        geo = data["Geo"]
        zAtm.dist = geo[0]
        zAtm.ang = geo[1]
        zAtm.tor = geo[2]
    def pushAtom(zAtm):
        zAtm.fixvarNum = len(atoms) + 4
        atoms.append(zAtm)
        zAtomDict[(resNum,atmName)] = zAtm
    cTerminal = (taggedSequence[len(taggedSequence)-1] == CTerminalChar)
    print(cTerminal)
    sequence = taggedSequence.replace(CTerminalChar,"")
    seqResidues = [None]
    for letter in sequence:
        seqResidues.append(PDBTools.Peptide.OneToThreeLetter[letter])
    print(seqResidues[1:])
    atoms = []
    zAtomDict = {}
    bbCountdown = 3 * centerResidue
    highestResidue = len(seqResidues) - 1
    threeQueue = [1,2,3]
    for resNum in range(centerResidue,0,-1):
        for atmName in SidechainData["reverseBackbone"]["Order"]:
            zAtom = popAtom(resNum,atmName)
            zAtom.backboneOrder = bbCountdown
            bbCountdown-=1
            applyGeo(zAtom,SidechainData["reverseBackbone"][atmName])
            pushAtom(zAtom)
            applyThreeQueue(zAtom,threeQueue)
            print(zAtom.toZMATLine())
    threeQueue = [4,5,6]
    for resNum in range(centerResidue+1,highestResidue+1):
        for atmName in SidechainData["forwardBackbone"]["Order"]:
            zAtom = popAtom(resNum,atmName)
            applyGeo(zAtom,SidechainData["forwardBackbone"][atmName])
            pushAtom(zAtom)
            applyThreeQueue(zAtom,threeQueue)
            print(zAtom.toZMATLine())
    for resNum in range(1,highestResidue+1):
        resType = seqResidues[resNum]
        resData = SidechainData[resType]
        for atmName in ("H","O"):
            zAtom = popAtom(resNum,atmName)
            if (atmName == "H"):
                if not resData["hasH"]:
                    continue
                zAtom.distAtom = zAtomDict[(resNum,"N")].fixvarNum
                zAtom.angAtom = zAtomDict[(resNum,"CA")].fixvarNum
                zAtom.torAtom = zAtomDict[(max(1,resNum-1),"C")].fixvarNum
                applyGeo(zAtom,SidechainData["H"])
                if (resNum == 1):
                    zAtom.tor =  0.00
            elif (atmName == "O"):
                zAtom.distAtom = zAtomDict[(resNum,"C")].fixvarNum
                zAtom.angAtom = zAtomDict[(resNum,"CA")].fixvarNum
                zAtom.torAtom = zAtomDict[(min(highestResidue,resNum+1),"N")].fixvarNum
                applyGeo(zAtom,SidechainData["O"])
                if (resNum == highestResidue):
                    zAtom.tor =  0.00
            pushAtom(zAtom)
    atmName = "CB"
    for resNum in range(1,highestResidue+1):
        resType = seqResidues[resNum]
        resData = SidechainData[resType]
        if not resData["hasCB"]:
            continue
        zAtom = popAtom(resNum,atmName)
        applyConnect(zAtom,SidechainData["CB"],resNum)
        pushAtom(zAtom)
    for resNum in range(1,highestResidue+1):
        resType = seqResidues[resNum]
        resData = SidechainData[resType]
        for atmName in resData["SideOrder"]:
            atmData = resData[atmName]
            zAtom = popAtom(resNum,atmName)
            applyConnect(zAtom,resData[atmName],resNum)
            pushAtom(zAtom)
    if cTerminal:
        zAtom = popAtom(highestResidue,"OXT")
        applyConnect(zAtom,SidechainData["OXT"],highestResidue)
        pushAtom(zAtom)
        zAtomDict[(highestResidue,"O")].torAtom = 1
    if outputFilename is None:
        outputFilename = "%s_aa%i.zmat"%(taggedSequence,centerResidue)
    makeZMATFile(atoms,outputFilename)
def makeZMATFile(atoms,filename):
    outf = open(filename,'w')
    outf.write("%6i MUL\n"%(len(atoms)))
    for a in atoms:
        print(a.toZMATLine())
        outf.write("%s\n"%a.toZMATLine())
    outf.close()
def applyThreeQueue(zAtm,threeQ):
    zAtm.distAtom = threeQ[0]
    zAtm.angAtom = threeQ[1]
    zAtm.torAtom = threeQ[2]
    threeQ[2] = threeQ[1]
    threeQ[1] = threeQ[0]
    threeQ[0] = zAtm.fixvarNum
    
if __name__ == "__main__":
    SeqToZMAT("GQTSV",3,"test.zmat")
'''
b = ZMAT.ZMAT_Atom.blank()
b.resType = "LYS"
print(b.identifier())
'''
'''
zm = ZMAT.ZMAT("TSV_aa3.zmat")
makeZMATFile(zm.zAtoms,"Test.zmat")
print(zm.zAtoms[4].identifier())
print(zm.zAtoms[4].toZMATLine())
print("     8 CA      7     6     5    1.51  120.00 -180.00 0 0 SER      2     5")
'''
