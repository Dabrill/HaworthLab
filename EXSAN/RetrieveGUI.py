import tkinter as tk
from tkinter import filedialog
import PDBTools
import Permutations
import os
targetProtein = None
targetLigand = None
def readTPFile(filepath,chains,xrayChain):
    global targetProtein
    global targetLigand
    inf = open(filepath,'r')
    tp = []
    lig = []
    for l in inf:
        a = PDBTools.Atom(l)
        if (a.valid):
            if (a.residueType not in ["HOH","WAT"]):
                if (a.chain in chains):
                    tp.append(a)
                elif  (xrayChain is None) or (a.chain in xrayChain):
                    a.chain = "X"
                    lig.append(a)
    targetProtein = PDBTools.Peptide(tp)
    targetLigand = PDBTools.Peptide(lig)
def getPoseHeaderNumber(line):
    if (len(line) == 0):
        return None
    if (line[0] == "["):
        if (line[-2:] == "]\n"):
            return int(line[1:-2])
    return None
def grabPoseScoreDataFromFile(poses,file):
    poseSet = set(poses)
    ret = [None] * len(poses)
    inf = open(file,'r')
    line = "initial"
    found = 0
    keepMode = False
    while (line != ""):
        line = inf.readline()
        poseNum = getPoseHeaderNumber(line)
        if (poseNum is None):
            if (keepMode):
                ret[i]+=line
        else:
            if poseNum in poseSet:
                keepMode = True
                found+=1
                i = poses.index(poseNum)
                ret[i] = ""
            else:
                keepMode = False
                if (found == len(poses)):
                    return ret
                    
    return ret
def grabScoringData(poses):
    align = grabPoseScoreDataFromFile(poses,"alignDump.dtf")
    print("Loaded align")
    side = grabPoseScoreDataFromFile(poses,"sideDump.dtf")
    print("Loaded side")
    hbond = grabPoseScoreDataFromFile(poses,"HBondTable.dtf")
    print("Loaded hbond")
    if (os.path.isfile("contable.tsv")):
        eTable = open("conTable.tsv").read().split("\n")
    else:
        eTable = open("ETable.tsv").read().split("\n")
    print("Loaded Energy Table")
    outf = open("RequestScoreData.dtf",'w')
    outf.write("%s\n"%eTable[0])
    for i in range(len(poses)):
        outf.write("%s\n"%eTable[poses[i]])

    
    tables = [(align,"Alignment:"),(side,"\nSide-Chain:"),(hbond,"\nHydrogen Bond:")]
    for i in range(len(poses)):
        outf.write("[%i]\n"%poses[i])
        for t in tables:
            outf.write("%s\n"%t[1])
            if (t[0][i] is None):
                outf.write("Missing\n")
            else:
                outf.write(t[0][i])
    outf.close()
class Application(tk.Frame):              
    def __init__(self, master=None):
        tk.Frame.__init__(self, master)   
        self.grid()                       
        self.createWidgets()
        self.useFile = False
    def doTheThing(self):
        #print(str(self.water.get())+","+str(self.hydrogen.get())+","+str(self.fileType.get())+","+str(self.allSubRand.get()))
        hydrogens = self.hydrogen.get()
        allSubRand = self.allSubRand.get()
        fileType = self.fileType.get()
        water = self.water.get()
        includeTarget = self.includeTarget.get()
        includeXRay = self.includeXRay.get()
        makeScoreFile = self.makeScoreFile.get()
        if (water == 1):
            try:
                poseWater = Permutations.WaterPermutations.readCWF("Water.cwf")
            except:
                water = 2
        if (hydrogens == 1):
            try:
                ligands = Permutations.PeptidePermutations.readCPF("AllH.cpf")
            except:
                ligands = Permutations.PeptidePermutations.readCPF("All.cpf")
        elif (hydrogens == 2):
            ligands = Permutations.PeptidePermutations.readCPF("All.cpf")
        if (allSubRand == 1):
            poseList = list(range(1,1+ligands.getLength()))
        elif (allSubRand == 2):
            poseList = self.getPoseList()
        elif (allSubRand == 3):
            import random
            poseList = random.sample(range(1,1+ligands.getLength()),int(self.textInstr.get()))
        if (makeScoreFile == 3):
            grabScoringData(poseList)
            self.quit()
            return
        if (fileType == 1):
            if (len(poseList) > 1):
                outfile = open("RequestedPoses.pdb",'w')
            else:
                outfile = open("Pose"+str(poseList[0])+".pdb",'w')
        if (includeTarget == 1):
            if (fileType == 1):
                outfile.write(str(targetProtein)+"TER\n")
        if (includeXRay == 1):
            if (fileType == 1):
                outfile.write(str(targetLigand)+"TER\n")   
        for p in poseList:
            poseStr = ""
            poseStr+=ligands.getPosePDB(p,False)
            if (water == 1):
                poseStr+=poseWater.grabPoseString(p,hydrogens==1)
            poseStr+="TER"+str(p)+"\n"
            if (fileType == 1):
                outfile.write(poseStr)
            else:
                outfile = open("Pose"+str(p)+".pdb",'w')
                if (includeTarget == 1):
                    outfile.write(str(targetProtein))
                if (includeXRay == 1):
                    outfile.write(str(targetLigand))
                outfile.write(poseStr)
                outfile.close()
        if (fileType == 1):
            outfile.close()
        if (makeScoreFile == 1):
            grabScoringData(poseList)
        self.quit()
    def getPoseList(me):
        if (me.useFile):
            precede,extenstion = os.path.splitext(me.filePath)
            if (extenstion == ".pos"):
                return PDBTools.readGoodPoseFile(me.filePath)
            else:
                contents = open(me.filePath,'r').read().split("\n")
                ret = []
                for c in contents:
                    try:
                        n = int(c)
                        ret.append(n)
                    except:
                        pass
                return ret
        else:
            split = me.textInstr.get().split(",")
            given = []
            for g in split:
                if (g.find("-") > -1):
                    halves = g.split("-")
                    for i in range(int(halves[0]),int(halves[1])+1):
                        given.append(i)
                else:
                    given.append(int(g))
            return given
    def createWidgets(self):
        self.goButton = tk.Button(self, text='Go',command=self.doTheThing)
        self.goButton.grid(row=0,column=0)   
        self.quitButton = tk.Button(self, text='Quit',command=self.quit)            
        self.quitButton.grid(row=0,column=1)
        self.water = tk.IntVar()
        self.water.set(2)
        #self.water = 1
        self.radios = []
        group = []
        group.append(tk.Radiobutton(self, text="Water", variable=self.water, value=1))
        group.append(tk.Radiobutton(self, text="No Water", variable=self.water, value=2))
        self.radios.append(group)
        self.hydrogen = tk.IntVar()
        self.hydrogen.set(2)
        group = []
        group.append(tk.Radiobutton(self, text="Hydrogen", variable=self.hydrogen, value=1))
        group.append(tk.Radiobutton(self, text="No Hydrogen", variable=self.hydrogen, value=2))
        self.radios.append(group)
        self.fileType = tk.IntVar()
        self.fileType.set(2)
        group = []
        group.append(tk.Radiobutton(self, text="Single File", variable=self.fileType, value=1))
        group.append(tk.Radiobutton(self, text="Separate Files", variable=self.fileType, value=2))
        self.radios.append(group)
        self.includeTarget = tk.IntVar()
        self.includeTarget.set(2)
        group = []
        group.append(tk.Radiobutton(self, text="Include Target", variable=self.includeTarget, value=1))
        group.append(tk.Radiobutton(self, text="Don't Include", variable=self.includeTarget, value=2))
        self.radios.append(group)
        self.includeXRay = tk.IntVar()
        self.includeXRay.set(2)
        group = []
        group.append(tk.Radiobutton(self, text="Include Xray Ligand", variable=self.includeXRay, value=1))
        group.append(tk.Radiobutton(self, text="Don't Include", variable=self.includeXRay, value=2))
        self.radios.append(group)
        self.allSubRand = tk.IntVar()
        self.allSubRand.set(2)
        group = []
        group.append(tk.Radiobutton(self, text="All Poses", variable=self.allSubRand, value=1))
        group.append(tk.Radiobutton(self, text="Subset", variable=self.allSubRand, value=2))
        group.append(tk.Radiobutton(self, text="Random", variable=self.allSubRand, value=3))
        self.radios.append(group)
        group = []
        self.makeScoreFile = tk.IntVar()
        self.makeScoreFile.set(2)
        group.append(tk.Radiobutton(self, text="Create Scoring File", variable=self.makeScoreFile, value=1, command = self.yesPDB))
        group.append(tk.Radiobutton(self, text="Don't Create", variable=self.makeScoreFile, value=2, command = self.yesPDB))
        group.append(tk.Radiobutton(self, text="Only (noPDB)", variable=self.makeScoreFile, value=3, command = self.noPDB))
        self.radios.append(group)


        numRows = len(self.radios)
        for row in range(numRows):
            group = self.radios[row]
            for i in range(len(group)):
                r = group[i]
                r.grid(row = 1 + row, column = i)
        self.textInstr = tk.StringVar()
        self.selections = tk.Entry(self,textvariable=self.textInstr, width=40)
        self.selections.grid(row=numRows+1,column=0,columnspan=3)
        self.fileDiag = tk.Button(self, text='Open File', command=self.askopenfilename)
        self.fileDiag.grid(row=numRows+2,column=0)
    def noPDB(me):
        for i in range(len(me.radios) - 2):
            group = me.radios[i]
            for button in group:
                button.config(state="disabled")
    def yesPDB(me):
        for i in range(len(me.radios) - 2):
            group = me.radios[i]
            for button in group:
                button.config(state="normal")
    def askopenfilename(me):
        path = os.path.abspath('')
        filename = filedialog.askopenfilename(initialdir = path,title = "Select file",filetypes = (("pose lists","*.pos"),("text file","*.txt"),("all files","*.*")))
        me.useFile = True
        me.filePath = filename
        dispName = filename
        index = dispName.find('/')
        while (index > -1):
            dispName = dispName[index+1:]
            index = dispName.find('/')
        me.textInstr.set(dispName)
        me.selections.config(state="disabled")
        me.allSubRand.set(2)
        for r in me.radios[4]:
            r.config(state="disabled")
        return filename
def getTPName():
    import json
    with open('../constants.json') as data_file:
        CONSTANTS = json.load(data_file)
        readTPFile(CONSTANTS["commonFolder"]+"/"+CONSTANTS["targetProtein"],CONSTANTS["targetProteinChain"],CONSTANTS.get("xrayChain",None))
TARGET_LOC = getTPName()
app = Application()                       
app.master.title('Retrieve TRUMP poses')    
app.mainloop()
