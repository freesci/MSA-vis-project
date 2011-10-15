import os
import fcntl

"""
Lock polegajacy na tym, ze dopoki plik djangolock.tmp jest w statusie locking (wywolane przez funkcje acquire)
 to inny proces nie moze odpalic na tym pliku funkcji acquire dopoki nie zostanie
wywolana funkcja release na pliku, odblokowujaca tenze plik
"""


#dokumentacja flock(2) -  http://fuse4bsd.creo.hu/localcgi/man-cgi.cgi?flock+2
class DjangoLock:

    def __init__(self, filename):
        self.filename = filename
        # This will create it if it does not exist already
        self.handle = open(filename, 'w')

    # flock() is a blocking call unless it is bitwise ORed with LOCK_NB to avoid blocking 
    # on lock acquisition.  This blocking is what I use to provide atomicity across forked
    # Django processes since native python locks and semaphores only work at the thread level
    def acquire(self):
	#fcntl.flock(handle, fcntl.LOCK_EX) jesli samo, to to wtedy process zablokowania jest zatrzymany, az nie zostanie rpez ten sam proces lub inny proces odblokowany
        fcntl.flock(self.handle, fcntl.LOCK_EX | fcntl.LOCK_NB)  #/* do not block when locking */
        #x | y   #Does a "bitwise or". Each bit of the output is 0 if the corresponding bit of x AND of y is 0, otherwise it's 1.

    def release(self):
        fcntl.flock(self.handle, fcntl.LOCK_UN)  #/* unlock file */

    def __del__(self):
        self.handle.close()


#USAGE:
#lock = DJangoLock(settings.MEDIA_ROOT + "uploaded_files/djangolock1.tmp")
#lock.acquire()
#try:
    #pass
#finally:
    #lock.release()
#del lock # == lock.__del__()