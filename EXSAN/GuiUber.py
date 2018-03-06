import tkinter as tk
from tkinter import messagebox
from tkinter import filedialog
from ctypes import c_char_p
from ctypes import c_bool
import os
import traceback
import time
import threading
import multiprocessing
import GU_Wrapper
import Log
LogFilename = "Megalog.txt"
root = os.path.abspath('')

        
class Application(tk.Frame):

    def startRun(self):
        Log.startTime = time.time()
        path = self.listbox.get(0)
        PRINT("\n\n\n\nNew Run on Inteface: %s"%path)
        
        
        self.good.value = True
        self.amRunning.value = True
        self.errorTrace.value = ""
        try:
            os.chdir(path)
            self.currentThread = multiprocessing.Process(target=GU_Wrapper.runExsan,args=(self.good,self.errorTrace))
            self.currentThread.start()
        except:
            self.good.value = False
            self.errorTrace.value = traceback.format_exc()
            self.currentThread = None
    def __init__(self, master=None):
        tk.Frame.__init__(self, master)
        self.grid()
        self.master.geometry('800x500')

        self.seconds = 0
        self.currentThread = None

        self.manager = multiprocessing.Manager()
        self.good = self.manager.Value(c_bool,False)
        self.amRunning = self.manager.Value(c_bool,False)
        self.errorTrace = self.manager.Value(c_char_p, "")
        
        self.createWidgets()
        self.useFile = False
    def createWidgets(self):
        self.goButton = tk.Button(self, text='Add',command=self.addRun)
        self.goButton.grid(row=0,column=0)   
        self.quitButton = tk.Button(self, text='Quit',command=self.quit)            
        self.quitButton.grid(row=0,column=1)

        self.timerText = tk.StringVar()
        self.timerText.set("0:00:00")
        self.timer = tk.Label(self.master, textvariable=self.timerText)
        self.timer.grid(row=1,column=0)
        
        self.listbox = tk.Listbox(self.master)
        self.listbox.config(width=130,height=20)
        self.listbox.grid(row=2,column=0)
        self.listbox.bind("<Double-Button-1>",self.removeLine)
        self.listbox.bind("<Up>",self.keyboardUpPress)
        self.listbox.bind("<Down>",self.keyboardDownPress)


    def keyboardUpPress(self,event):
        index = self.listbox.curselection()[0]
        text =  self.listbox.get(index)
        if (index > 1):
            self.listbox.selection_clear(index)
            self.listbox.delete(index)
            self.listbox.insert(index-1, text)
            self.listbox.activate(index)
    def keyboardDownPress(self,event):
        index = self.listbox.curselection()[0]
        text =  self.listbox.get(index)
        if (index > 0):
            self.listbox.selection_clear(index)
            self.listbox.delete(index)
            self.listbox.insert(index+1, text)
            self.listbox.activate(index)        
    def removeLine(self,event):
        index = self.listbox.curselection()[0]
        selection = self.listbox.get(index)
        if (index > 0):
            confirmed =  messagebox.askokcancel("Confirm","Remove Run %s?"%selection)
            if confirmed:
                self.listbox.delete(index)
        else:
            confirmed =  messagebox.askokcancel("Halt Run While In Process?","Halt Run %s?"%selection)
            if confirmed:
                self.listbox.delete(index)
                self.abort()
    def nextRun(self):
        self.listbox.delete(0)
        self.seconds = 0
        if (self.listbox.size() == 0):
            PRINT("Queue Offline")
            self.amRunning.value = False
        else:
            self.startRun()
    def update(self):
        if not self.amRunning.value:
            if (self.listbox.size() > 0):
                PRINT("Queue Online")
                self.startRun()
        if self.amRunning.value:
            self.seconds+=1
            if (self.currentThread is None) or (not self.currentThread.is_alive()):
                #self.currentThread.join()
                PRINT("Complete")
                if self.good.value:
                    PRINT("Success")
                    PRINT(secondsToHMS(self.seconds))
                else:
                    PRINT("Run failed")
                    PRINT(self.errorTrace.value)
                os.chdir(root)
                self.nextRun()
        self.timerText.set(secondsToHMS(self.seconds))
        self.master.after(1000,self.update)
    def abort(self):
        self.currentThread.terminate()
        self.good.value = False
        self.errorTrace.value = "Run Halted By User!"
        self.amRunning.value = False        
    def addRun(self):
        path = os.path.abspath('')
        startPath = root
        startPath = "D:\\Exsan\\PDZ"
        filename = tk.filedialog.askdirectory(initialdir = startPath,title = "Select folder for run")
        if isARunFolder(filename):
            self.listbox.insert(tk.END, filename)
        else:
            messagebox.showerror("Not a run folder","Error: Selected folder is not a run folder")
    def quit(self):
        if self.amRunning.value:
            confirmed =  messagebox.askokcancel("Quit and end run?","Are you want to quit while run is active?")
            if confirmed:
                self.abort()
            else:
                return
        self.master.destroy()
def isARunFolder(path):
    return os.path.isfile("%s/plan.json"%path) and os.path.isfile("%s/constants.json"%path)
def secondsToHMS(secs):
    return "%i:%02i:%02i"%(secs/3600,(secs%3600)/60,secs%60)
def PRINT(s):
    logFile = open("%s/%s"%(root,LogFilename),'a')
    print(s)
    logFile.write(s+"\n")
    logFile.close()
def Main():
    PRINT("Booting Queue Gui\n\n")
    app = Application()
    '''
    app.listbox.insert(tk.END,"D:/Exsan/PDZ/2QBW/Xray")
    app.listbox.insert(tk.END,"D:/Exsan/PDZ/BE12/Xray")
    app.listbox.insert(tk.END,"D:/Exsan/PDZ/4Z8J/Xray")
    app.listbox.insert(tk.END,"D:/Exsan/PDZ/5I7Z/Xray2")
    app.listbox.insert(tk.END,"D:/Exsan/PDZ/BE9/Testbed")
    '''
    app.master.title('EXSAN Queue Manager')
    app.master.after(10,app.update)
    app.mainloop()
if __name__ == "__main__":
    Main()



