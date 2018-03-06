import Vector
import math
import copy
import os
def invertDictionary(dictIn):
    dictOut = {}
    count = 0
    for entry in dictIn:
        key = dictIn[entry]
        if (key not in dictOut):
            dictOut[key] = entry
            count+=1
        else:
            raise Exception("Input dictionary is not one to one")
    return dictOut
class Peptide:
    ThreeToOneLetter = {"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLU":"E","GLN":"Q","GLY":"G","HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V"}
    OneToThreeLetter = invertDictionary(ThreeToOneLetter)
    def __init__(self,atoms,indexNum = 0,CAList = None):
        self.atoms = atoms
        self.indexNum = indexNum
        if (CAList == None):
            self.indexCA()
        else:
            self.CAList = CAList
        self.calcSequence()
        self.residues = None
        self.firstRes = None
    @staticmethod
    def fromStringArray(array,index = None):
        atoms = []
        for l in array:
            if "TER" not in l:
                atoms.append(Atom(l))
            else:
                if (len(l) > 3):
                    index = int(l[3:])
        return Peptide(atoms,index)
    def indexCA(self):
        self.CAList = []
        for a in self.atoms:
            if (a.atomType == "CA"):
                self.CAList.append(a)
    def listResidues(me):
        residues = []
        lastAdd = None
        for a in me.atoms:
            if (lastAdd is None) or (a.residueNumber != lastAdd):
                residues.append((a.residueNumber,a.residueType))
                lastAdd = a.residueNumber
        return residues
    def calcSequence(me):
        residues = me.listResidues()
        me.sequence = ""
        for r in residues:
            me.sequence+=Peptide.ThreeToOneLetter.get(r[1],"X")        
    def makeDictionary(me):
        me.dictOffset = me.atoms[0].residueNumber
        me.highestResidue = max(me.atoms,key=lambda x:x.residueNumber).residueNumber
        me.residues = []
        for i in range(me.highestResidue - me.dictOffset + 1):
            me.residues.append({})
        for a in me.atoms:
            me.residues[a.residueNumber - me.dictOffset][a.atomType] = a
    def getResidue(me,residueNum):
        return me.residues[residueNum - me.dictOffset]                
    def getSpecificAtom(me,residueNum,atmType):
        res = me.residues[residueNum - me.dictOffset]
        if atmType in res:
            return res[atmType]
        return None
    def getAtomNumber(me,number):
        if (number < len(me.atoms)):
            if (me.atoms[number-1].atomNumber == number):
                return me.atoms[number-1]
        for a in me.atoms:
            if (a.atomNumber == number):
                return a
        return None
    def getSequence(me):
        return me.sequence
    def findOverlaps(big,small):
        if isinstance(small,str):
            smallSeq = small
        elif isinstance(small,Peptide):
            smallSeq = small.getSequence()
        else:
            raise Exception("Must input string or peptide object")
        return big.listAllSubstringHits(big.sequence,smallSeq)
            
    @staticmethod
    def listAllSubstringHits(base,key):
        pos = 0
        hits = []
        while (pos > -1):
            pos = base.find(key,pos)
            if (pos > -1):
                hits.append(pos)
                pos+=1
        return hits        
    def toTemplate(me):
        newAtoms = []
        for i in range(len(me.atoms)):
            a = copy.copy(me.atoms[i])
            coords = []
            for d in range(3):
                coords.append(a.location.dims[d]*1000.0)
            a.location = copy.copy(a.location)
            a.location.dims = coords
            newAtoms.append(a)
        return Peptide(newAtoms,me.indexNum)
    def fromTemplate(me):
        newAtoms = []
        for i in range(len(me.atoms)):
            a = copy.copy(me.atoms[i])
            coords = []
            for d in range(3):
                coords.append(a.location.dims[d]/1000.0)
            a.location = copy.copy(a.location)
            a.location.dims = coords
            newAtoms.append(a)
        return Peptide(newAtoms,me.indexNum)
    def minDistance(self,other):
        minCA = 9999999
        for i in range(0,len(self.CAList)):
            a = self.CAList[i]
            resMin = 999999
            for j in range(0,len(other.CAList)):   
                b = other.CAList[j]
                dist = a.distance(b)
                if (dist < minCA):
                    minCA = dist           
        return minCA
    def minDistanceAtom(self,other):
        minAtm = self.atoms[0].distance(other)
        for i in range(0,len(self.atoms)):
            a = self.atoms[i]
            dist = a.distance(other)
            if (dist < minAtm):
                minAtm = dist           
        return minAtm
    def sameAtom(self,other,atom):
        for i in range(0,3):
            if (not (self.atoms[atom][i] - other.atoms[atom][i] == 0)):
                return False
        return True
    def resInDistance(self,tree,distance,residue):
        best = None
        for a in self.atoms:
            if (a.residueNumber == residue):
                closest, dist = tree.nearestNeighbor(a.location,True)
                if best is None:
                    best = dist
                else:
                    if (dist < best):
                        best = dist
        if best is None:
            return False
        if (best > distance):
            return False
        return True
    def allowedDistance(self,tree,minNeeded,maxAtResidue,numException,allAtoms = False):
        if (allAtoms):
            useList = self.atoms
        else:
            useList = self.CAList
        minFound = False
        numberOfAtoms = len(self.atoms)
        highRes = self.atoms[numberOfAtoms-1].residueNumber
        lowRes = self.atoms[0].residueNumber
        numOfRes = highRes - lowRes + 1
        #print("numberOfResidues ",numOfRes)
        bests = [None] * numOfRes

        
        for i in range(0,len(useList)):
            a = useList[i]
            res = a.residueNumber
            closest = tree.nearestNeighbor(a.location)
            dist = closest.distance(a)
            if (dist <= minNeeded):
                minFound = True
            if bests[res-lowRes] is None:
                bests[res-lowRes] = dist
            else:
                if (dist < bests[res-lowRes]):
                    #print(a.identifier(),closest.identifier())
                    bests[res-lowRes] = dist

        exception = 0
        if not minFound:
            return False
        #print("Bests: ",bests)
        for b in bests:
            if (b > maxAtResidue):
                exception+=1
        if (exception > numException):
            return False
        return True  
    def toMatrix(self):
        ret = []
        for a in self.atoms:
            ret.append(a.location.dims)
        return ret
    def toIntMatrix(self):
        ret = []
        for a in self.atoms:
            e = [None,None,None]
            for d in range(3):
               #print(a.location.dims[d])
               s = str(a.location.dims[d])
               #print(s)
               e[d] = int(s.replace(".",""))
               #print(e[d])
            ret.append(e)
        return ret
    def simplifySelf(self):
        CAList = None
        for i in range(0,len(self.atoms)):
            self.atoms[i] = self.atoms[i].location.dims
    def __str__(self):
        return self.makePDBString(False)
    def makePDBString(self,addTer = True):
        s = ""
        for a in self.atoms:
            s+=a.toPDBLine()+"\n"
        if (addTer):
            s+="TER"+str(self.indexNum)+"\n"
        return s
    def toStringArray(self,addTer = True):
        s = []
        for a in self.atoms:
            s.append(a.toPDBLine())
        return s
    def compareMissing(self,other,start):
        resIndx = -1
        lastRes = None
        total = 0
        count = 0.0
        missing = []
        for a in self.atoms:
            if (lastRes is None) or (a.residueNumber > lastRes):
                resIndx+=1
                lastRes = a.residueNumber
            match = None
            try:
                match = other.residues[resIndx+start][a.atomType]
            except:
                missing.append(a)
        return missing
    def compare(self,other,start = None):
        if other.residues is None:
            other.makeDictionary()
        if start is None:
            overlaps = other.findOverlaps(self)
            if (len(overlaps) == 0):
                raise Exception("Sequences do not overlap")
            start = overlaps[0]
        resIndx = -1
        lastRes = None
        total = 0.0
        backbone = 0.0
        count = 0.0
        bbCount = 0.0
        for a in self.atoms:
            if (lastRes is None) or (a.residueNumber > lastRes):
                resIndx+=1
                lastRes = a.residueNumber
            match = None
            try:
                match = other.residues[resIndx+start][a.atomType]
                #print(match.identifier())
                dist = a.distance(match)**2
                total+=dist
                count+=1
                if (a.isBackbone()):
                    backbone+=dist
                    bbCount+=1
            except:
                pass
        return total / count, backbone / bbCount
    def copy(me):
        import copy
        ret = Peptide([])
        ret.indexNum  = me.indexNum
        ret.residues = []
        for a in me.atoms:
            ret.atoms.append(copy.copy(a))
        return ret
    def createNoHCopy(me):
        raise NotImplementedError
        ret = Peptide([])
        ret.indexNum  = me.indexNum
        count = 1
        for a in me.atoms:
            if (a.isHydrogen()):
                print("\t"+str(a))
            else:
                print(str(a))
    def getCTerminalCarbon(me):
        last = None
        hasOxt = False
        for a in reversed(me.atoms):
            if (a.atomType == "OXT"):
                hasOxt = True
            if (a.atomType == "C") and (hasOxt):
                a.cTerminalCarbon = True
                return a
            if (last is None):
                last = a.residueNumber
            if (a.residueNumber < last) and (not a.isHydrogen()):
                return None
        return None
    def implyHydrogen(me,res):
        try:
            if (me.residues is None):
                me.makeDictionary()
            h = me.getSpecificAtom(res,"H")
            if (h is not None):
                return
            h = me.getSpecificAtom(res,"HN")
            if (h is not None):
                return h
        except:
            print("Crash ",res)
        if (res - me.dictOffset < 0):
            raise Exception("Residue number too low in imply Hydrogen: %i/%i"%(res,me.dictOffset))
        if (res - me.dictOffset >= len(me.residues)):
            raise Exception("Residue number too high in imply Hydrogen: %i/%i"%(res,me.dictOffset+len(me.residues)))       
        CA = me.getSpecificAtom(res,"CA")
        N = me.getSpecificAtom(res,"N")
        if (res - me.dictOffset > 0):
           C = me.getSpecificAtom(res-1,"C")
           H = Atom.getH(N,CA,C)
        else:
           C = me.getSpecificAtom(res,"C")
           H = Atom.getH(N,CA,C,True)
        me.residues[res - me.dictOffset]["H"] = H
        me.atoms.append(H)
        return H
    def getTerminalOxygen(pep):
        if (pep.residues is None):
            for a in pep.atoms:
                if (a.atomType == 'OXT'):
                    return a
            return None
        lastResidue = pep.atoms[len(pep.atoms)-1].residueNumber
        return pep.getSpecificAtom(lastResidue,"OXT")

            
            
        
class Atom:
    hydrophilicTable = {
             "ARG":"CZ",
             "LYS":"NZ",
             "ASP":"CG",
             "GLU":"CD",
             "GLN":"CD",
             "ASN":"CG",
             "HIS":"CE1",
             "SER":"OG",
             "THR":"OG1",
             "TYR":"OH"}
    hydrophobicTable = {
             "ALA":("CB",),
             "VAL":("CB","CG1","CG2"),
             "LEU":("CB","CG","CD1","CD2"),
             "ILE":("CB","CG1","CG2","CD"),
             "MET":("CB","CG","SD","CE"),
             "CYS":("CB","SG"),
             "PHE":("CB","CG","CD1","CD2","CE1","CE2","CZ"),
             "PRO":("CB","CG","CD"),
             "TYR":("CB","CG","CD1","CD2","CE1","CE2","CZ"),
             "TRP":("CB","CD2","CE2","CE3","CZ2","CZ3","CH2","CD1","CG")}
    ABtable = {"GLU":"A",
             "ASP":"A",
             "LYS":"B",
             "ARG":"B"}
    backboneTable = ("C","CA","N","OXT")
    def __init__(self,ln):
        ltype = ln[:6]
        isHetero = (ltype == "HETATM") or (ltype == "HETA  ")
        if (ltype == "ATOM  ") or isHetero:
            self.atomNumber = int(ln[7:11])
            self.atomType = ln[12:16].replace(" ","")
            self.residueType = ln[17:20].replace(" ","")
            self.chain = ln[21:22].replace(" ","")
            self.residueNumber = int(ln[22:26])
            self.location = Vector.Vector3D(float(ln[31:38]),float(ln[38:46]),float(ln[46:54]))
            self.valid = True
            self.hetero = isHetero
        else:
            self.valid = False
            #print (ln)
    @staticmethod
    def blankAtom():
        a = Atom("")
        a.valid = False
        return a
    def copy(self,copyLoc = True):
        a = Atom("")
        a.atomNumber = self.atomNumber
        a.atomType = self.atomType
        a.residueType = self.residueType
        a.chain = self.chain
        a.residueNumber = self.residueNumber
        a.valid = self.valid
        a.hetero = getattr(self, 'hetero',False)
        if (copyLoc):
            a.location = self.location.copy()
        else:
            a.location = None
        return a
            
        
    @staticmethod
    def getH(N,CA,C,internal = False):
        H = Atom.blankAtom()
        H.atomNumber = 0
        H.atomType = "H"
        H.residueType = N.residueType
        H.residueNumber = N.residueNumber
        H.chain = N.chain
        if internal:
            H.setAtomLocationByZMAT(N,CA,C,1.01,120,0)
        else:
            H.setAtomLocationByZMAT(N,C,CA,1.01,120,180)
        return H
    def setAtomLocationByZMAT(me,aDist,aAng,aTor,dist,ang,tor):
        #print(dist,ang,tor)
        me.location = Vector.Vector3D.zmatToCartesian(aDist.location,aAng.location,aTor.location,dist,ang,tor)
        for i in range(3):
            me.location.dims[i] = round(me.location.dims[i],3)
    def isCTerminalCarbon(atm):
        try:
            return atm.cTerminalCarbon
        except:
            return False
    def acidBase(atm):
        if atm.isCTerminalCarbon():
            return "A"
        if atm.residueType in Atom.ABtable:
            return Atom.ABtable[atm.residueType]
        return "N"
    def isBackbone(atm):
        if atm.atomType in Atom.backboneTable:
            return True
        return False
    def getHydrophilicAtom(atm):
        return Atom.hydrophilicTable.get(atm.residueType,None)
    def isHydrophilicRes(atm):
        if atm.residueType in Atom.hydrophilicTable:
            return True
        return False
    def isHydrophilicAtom(atm):
        if (atm.isCTerminalCarbon()):
            return True
        keyAtom = Atom.hydrophilicTable.get(atm.residueType,None)
        return (keyAtom == atm.atomType)
    def isHydrophobicRes(atm):
        if atm.residueType in Atom.hydrophobicTable:
            return True
        return False
    def isHydrophobicAtom(atm):
        aminoAcid = Atom.hydrophobicTable.get(atm.residueType,None)
        if aminoAcid is None:
            return False
        if atm.atomType in aminoAcid:
            return True
        return False
    def isHydrogen(self):
        if (len(self.atomType) > 3):
            return True
        if (self.atomType[0] == "H"):
            return True
        return False
    def distance(self,other):
        return self.location.distance(other.location)
    def angle(self,dist,ang,roundAnswer=True):
        ang = self.location.angleThreePoints(dist.location,ang.location)
        if (roundAnswer):
            return Vector.cleanAngle(ang)
        return ang
    def torsion(self,dist,ang,tor,roundAnswer = True):
        return self.location.dihedral(dist.location,ang.location,tor.location,roundAnswer)
    def toPDBLine(a):
        if (len(a.atomType)  < 4):
            aType = " "+a.atomType
        else:
            aType = a.atomType
        if getattr(a, 'hetero',False):
            lineType = "HETATM"
        else:
            lineType = "ATOM"
        return "%-6s%5i %-4s %s %s%4i    % 8.3f% 8.3f% 8.3f"%(lineType,a.atomNumber,aType,a.residueType,a.chain,a.residueNumber,a.location.dims[0],a.location.dims[1],a.location.dims[2])
    def resStr(a):
        if (a.residueType in ["WAT","HOH"]):
            return "%s%d"%(a.residueType,a.residueNumber)
        return "%s-%s%d"%(a.chain,a.residueType,a.residueNumber)
    def identifier(a):
        if (a.residueType in ["WAT","HOH"]):
            return "%s-%d %s"%(a.residueType,a.residueNumber,a.atomType)
        return "%s-%s-%d %s"%(a.chain,a.residueType,a.residueNumber,a.atomType)
        #return a.residueType+"-"+str(a.residueNumber)+" "+a.atomType
    def __str__(self):
        return self.identifier()
def toPDBFormat(vector):
    return formatNumber(vector.dims[0])+" "+formatNumber(vector.dims[1])+" "+formatNumber(vector.dims[2])
def readAtomsInChain(file,chain,aType = None,returnPeptide = False):
    #print(file)
    import time
    inFile = open(file, 'r')
    table=[]
    for line in inFile:
        line = line.replace('\n','')
        if(line[:4] == "ATOM") or (line[:6] == "HETATM"):
            atm = Atom(line)
            if (atm.chain == chain):
                if aType is not None:
                    if (atm.atomType[:len(aType)] == aType):
                        table.append(atm)
                else:
                    table.append(atm)
    if (returnPeptide):
        return Peptide(table,1)
    return table
def readPDBSelectively(file,chain = None,aType = None,includeWaters = False,includeHetero = False):
    if (aType is not None) and (not isinstance(aType,(list,tuple))):
        aType = (aType,)
    inFile = open(file, 'r')
    table=[]
    #print("Chain: ",chain)
    for line in inFile:
        line = line.replace('\n','')
        atm = Atom(line)
        if atm.valid:
            if (not includeWaters) and( atm.residueType in ["HOH","WAT"]):
                continue
            if (chain is not None) and (atm.chain not in chain):
                continue
            if (aType is not None) and (atm.atomType not in aType):
                continue
            if (not includeHetero) and (atm.hetero):
                continue
            table.append(atm)
    return table
def readPDBFile(file,excludeWater = True):
    inFile = open(file, 'r')
    table=[]
    newMol = []
    mol = 1
    CAList = []
    for line in inFile:
        line = line.replace('\n','')
        if ("ATOM" == line[:4]):
            atm = Atom(line)
            good = True
            if (excludeWater):
                good = good and (not atm.residueType in ["HOH","WAT"])
            if (good):
                newMol.append(atm)
            if (atm.atomType == "CA"):
                CAList.append(atm)
        if ("TER" == line[:3]):
            if (len(newMol) > 0):
                if (len(line) > 3):
                    try:
                        fileNum = int(line[3:])
                    except:
                        fileNum = None
                    table.append(Peptide(newMol,fileNum,CAList))
                else:
                    table.append(Peptide(newMol,mol,CAList))
                    mol+=1
                newMol = []
                CAList = []
    inFile.close()
    if (len(newMol) > 0):
        table.append(Peptide(newMol,mol,CAList))
    return table
def readAllInPDB(file):
    inf = open(file,'r')
    ret = [Atom(line) for line in inf.readlines()]
    inf.close()
    return ret
def fileToStringArray(file):
    fileRead = open(file,'r')
    back = []
    for line in fileRead:
        back.append(line)
    fileRead.close()
    return back   
def StringArrayToString(array):
    s = ""
    for line in array:
        s+=line
        if not line[-1:] == "\n":
            #print("|"+line[-1:]+"|")
            s+="\n"          
    return s
def newLoc(angle):
    dx = math.cos(angle)
    dy = math.sin(angle)
def padString(s,numChars,front = True):
    r = s
    while (len(r) < numChars):
        if (front):
            r = " "+r
        else:
            r = r+" "
    return r
def floatToInt(d):
    s = formatNumber(d)
    return(int(s.replace(".","")))
def formatNumber(x):
    s = '%07.3f' % (x)
    return s
def foundStreak(inAry,start):
    cur = start
    while (len(inAry) - 1 > cur) and (inAry[cur+1] == inAry[cur] + 1):
        cur = cur + 1
    return (cur - start)
def readGoodPoseFile(filename):
    outarray = []
    inf = open(filename,'r')
    pastHeader = False
    for line in inf:
        if (pastHeader):
            instruction = line[:1]
            if (instruction == "#"):
                outarray.append(int(line[1:]))
            else:
                parameters = line[1:].split(",")
                parameters = [int(x) for x in parameters]
                if (instruction == "+"):
                    for n in range(parameters[0],parameters[0]+parameters[1]+1):
                        outarray.append(n)
                if (instruction == "-"):
                    for n in range(parameters[0],parameters[1]+1):
                        outarray.append(n)
        else:
            pastHeader = True
    return outarray
def dumpArray(file,array):
    outf = open(file,'w')
    for n in array:
        outf.write(str(n)+"\n")
    outf.close()
def writeGoodPoseFile(filename,inAry):
    outf = open(filename,'w')
    cur = 0
    CUTOFF = 10
    outf.write(str(len(inAry))+"\n")
    while (cur < len(inAry)):
        streak = foundStreak(inAry,cur)
        if (streak == 0):
            outf.write("#"+str(inAry[cur])+"\n")
            cur+=1
        elif (streak < CUTOFF):
            outf.write("+"+str(inAry[cur])+","+str(streak)+"\n")
            cur+=streak+1
        elif (streak >= CUTOFF):
            outf.write("-"+str(inAry[cur])+","+str(inAry[cur]+streak)+"\n")
            cur+=streak+1
        else:
            raise Exception("Serious Error. Streak value not valid: "+str(streak))
def pickOverlap(overlaps,anchor):
    if (len(overlaps) == 0):
        return -1
    elif (len(overlaps) == 1):
        if anchor is None:
            return overlaps[0]
        if anchor.upper() in ["C","N"]:
            return overlaps[0]
        if (anchor == overlaps[0]):
            return overlaps[0]
        return -3    
    elif (len(overlaps) > 1):
        if isinstance(anchor,str):
            if (anchor.upper() ==  "C"):
                return overlaps[len(overlaps)-1]
            elif (anchor.upper() == "N"):
                return overlaps[0]
        elif isinstance(anchor,int):
            if anchor in overlaps:
                return anchor
        return -2
    return -4  
def lastResidueNum():
    currentFolder = os.path.abspath('')
    base = readPDBFile("out.pdb")
    return base[0].atoms[len(base[0].atoms)-1].residueNumber
def commaList(myList):
    ret = myList[0]
    for i in range(1,len(myList)):
        ret+=(',%s'%myList[i])
    return ret
def readHetWaters(filename):
    inf = open(filename,'r')
    atoms = []
    for l in inf:
        if (l[:6] == "HETATM"):
            if (l[17:20] == "HOH"):
                a = Atom(l)
                atoms.append(a)
    return atoms

