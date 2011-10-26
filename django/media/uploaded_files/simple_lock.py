import os,fcntl

#dokumentacja flock(2) -  http://fuse4bsd.creo.hu/localcgi/man-cgi.cgi?flock+2
class DjangoLock:

    def __init__(self, filename):
        self.filename = filename
        # This will create it if it does not exist already
        self.handle = open(filename, 'w')

    def acquire(self):
        fcntl.flock(self.handle, fcntl.LOCK_EX)  #/* exclusive file lock */

    def release(self):
        fcntl.flock(self.handle, fcntl.LOCK_UN)  #/* unlock file */

    def __del__(self):
        self.handle.close()