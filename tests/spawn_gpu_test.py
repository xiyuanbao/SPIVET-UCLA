import multiprocessing
def spawn(x):
    print 'Spawned',x
    import openpiv.gpu_process
    
    runpiv()
    
def runpiv():
    
    import openpiv.gpu_process
    print "here"
class test():
    def run(self):
        procs=[]
        for i in range(5):
            x=i+1
            p = multiprocessing.Process(target=spawn,args=(x,))
            p.start()
            procs.append(p)
        [p.join() for p in procs]
t=test()
t.run()
