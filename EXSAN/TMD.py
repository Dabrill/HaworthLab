import os
import shutil
def runTMD(CONSTANTS,step,fxvr = None,recorder = None):
    def PRINT(string):
        print(string)
        if recorder is not None:
            recorder.write(string+"\n")
    PRINT("Starting TMD run")
    poses = 0
    currentFolder = os.path.abspath('')
    commonFolder = CONSTANTS["commonFolder"]
    targetFilename = CONSTANTS["targetProtein"]
    if (fxvr is not None):
        fxvr.createApprovedFixvar("fixvar.in")
        #shutil.copy(dname+"/Fixvar.out","fixvar.in")
        posesIn = fxvr.numGood()
            
        #posesIn = poseNumRecorder.lastNumberPoses()
        PRINT(str(posesIn)+" poses input")
    else:
        PRINT("Starting virgin run")
        posesIn = "Nonsense"
    createDataDotIn(CONSTANTS,step,posesIn)
    
    #Checking to see if peptide is adding on the reverse end
    #Since it affects making the PDB file
    
    reverse = step["reverse"]
    
    shutil.copy(commonFolder+"/"+targetFilename,currentFolder)
    for file in os.listdir(commonFolder):
        if file.endswith(".zmat"):
            shutil.copy(commonFolder+"/"+file,currentFolder)
    PRINT("Starting TMD Calculations")

    os.system("%s/%s"%(CONSTANTS["universalFolder"],CONSTANTS["TMDLoc"]))
    PRINT("TMD program complete")
    summary = open("data.add","r")
    for line in summary:
        if "number of solutions:" in line:
            i = line.index(':')
            x = line[i+1:]
            poses = int(x)
                

    
    #PRINT("Making directories")

    try:
        os.remove("fixvr.in")
    except:
        pass
    try:
        os.remove("go_dock.bat")
        os.remove("go_intra.bat")
        os.remove("go_makeonepdb.bat")
        os.remove("go_makepdb.bat")
        os.remove("intra.in")
        os.remove("dock.in")
        os.remove(targetFilename)
    except:
        PRINT("No TMD ouput. Zero solution poses probable.")
        
    if (poses > 0):
        os.remove("makepdb.in")
        #os.remove("makeonepdb.in")

    
    shutil.copy("data.in","data.old")
    os.remove("data.in")
    removeWithKey(currentFolder,".zmat")
    PRINT("TMD Done")
    return poses
def referenceAtom(CONSTANTS):
    return padString(str(CONSTANTS["referenceAtomDistance"]),5)+padString(str(CONSTANTS["referenceAtomAngle"]),5)+padString(str(CONSTANTS["referenceAtomTorsion"]),5)
def createDataDotIn(CONSTANTS,step,posesIn):
    targetFilename = CONSTANTS["targetProtein"]
    outf = open("data.in",'w')
    outf.write(step["sequence"]+" to "+CONSTANTS["tpName"]+"\n")
    outf.write("  0  0  1  1  0  0  0  0\n")
    outf.write(padString(targetFilename,21,False)+CONSTANTS["targetProteinOptions"]+" "+CONSTANTS["targetProteinChain"]+"\n")
    outf.write(referenceAtom(CONSTANTS)+"\n")
    outf.write(step["zsequence"]+step["zmatSuffix"]+".zmat\n")
    outf.write("  1  0  0  0  0  0        idok,ilgout,nwropt,ntargp,iintml,isamin\n")
    outf.write(padString(str(len(step["variableTorsions"])),3)+padString(str(step["nclsrs"]),3)+padString(str(step["numFixvar"]),3)+padString(str(len(step["criteria"])),3)+padString(str(step["nDistanceCheck"]),3)+"\n")
    for i in range(0,len(step["variableTorsions"])):
        twist = step["variableTorsions"][i]
        s=padString(str(twist["atom"]),3)
        s+=padString(str(twist["start"]),4)+padString(str(twist["stop"]),4)+padString(str(twist["step"]),4)
        s+=padString(str(len(twist["HoldAngle"])),2)+padString(str(twist["type"]),2)+"\n"
        outf.write(s)
        for j in range(0,len(twist["HoldAngle"])):
            outf.write(padString(str(twist["HoldAngle"][j]["atom"]),3)+padString(str(twist["HoldAngle"][j]["angle"]),4)+"\n")      
    try:
        clashLimit = step["specialClash"]
    except:
        clashLimit = CONSTANTS["clashLimit"]
    outf.write("  %.1f   %i\n"%(clashLimit,CONSTANTS["allowedClash"]))
    #outf.write("  "+str(clashLimit)+"   "+str(CONSTANTS["allowedClash"])+"\n")
    if (step["numFixvar"]> 0):
        outf.write(padString(str(posesIn),8)+"    "+str(step["reverse"])+"\n")
    for i in range(0,len(step["criteria"])):
        crite = step["criteria"][i]
        criteS = "%5i %-3s  %3i %-3s%5i%8i%11.1f%7.1f\n"%(crite["refPNum"],crite["refPAtm"],crite["refTNum"],crite["refTAtm"],crite["param1"],crite["param2"],crite["minDist"],crite["maxDist"])
        outf.write(criteS)
        #outf.write(padString(str(crite["refPNum"]),5)+" "+padString(crite["refPAtm"],2,False)+padString(str(crite["refTNum"]),6)+" "+padString(crite["refTAtm"],2,False)+padString(str(crite["param1"]),6)+padString(str(crite["param1"]),8)+padString(str(crite["minDist"]),11)+padString(str(crite["maxDist"]),7)+"\n")
    outf.close()
def removeWithKey(dirName,key):
    for file in os.listdir(dirName):
        if file.endswith(key):
            os.remove(dirName+"/"+file)
def padString(s,numChars,front = True):
    r = s
    while (len(r) < numChars):
        if (front):
            r = " "+r
        else:
            r = r+" "
    return r
def cleanDirectories():
    removeFolder("makePdb")
    removeFolder("order")
    removeFolder("makeOnePdb")
    removeFolder("watgenorder")
def removeFolder(top):
    for root, dirs, files in os.walk(top, topdown=False):
        for name in files:
            os.remove(os.path.join(root, name))
        for name in dirs:
            os.rmdir(os.path.join(root, name))
