import traceback
from EXSAN_11_27_17 import Main
def runExsan(status,trace):    
    try:
        back = Main()
        if (back!=""):
            status.value = False
            trace.value = back
    except Exception:
        status.value = False
        trace.value = traceback.format_exc()
