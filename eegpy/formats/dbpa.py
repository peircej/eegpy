#! /usr/bin/python
# -*- coding: utf8 -*-
import os
import tempfile
import numpy as np

from eegpy.formats.iobase import (MemoryMappedBinaryDataFile,
                                  EEGfiltered)

class DBPA(MemoryMappedBinaryDataFile):
    """Access to Sensorium DBPA datafiles"""
    _eventChannels = [0]
    _responseChannels = [1] 
    _dataChannels =  range(2,119)
    _eyeChannels = range(119,122)

    def __init__(self,filename,mode="r+",cNames=None,shape=None, Fs = 1000.0, numChannels=122):
        """Documentation needed"""
        fileSize = os.stat(filename).st_size
        numDatapoints = int(fileSize/numChannels/4)
        self.numChannels = numChannels
        MemoryMappedBinaryDataFile.__init__(self, filename, mode, 
                    shape = (numDatapoints, numChannels),
                    data_offset=0,
                    dtype=">f", Fs = Fs, channel_names=cNames)
    def getEventIndices(self, onsets=True, eventChannel=0, excludeEvents=[0]):
        """Retrieve the indices of all the events and the array of what the IDs were
        starting/ending
        
        :params:
                onsets: True, False indicates whether to return onset events or offset (onsets=False)
                eventChannel: int to specify which of the event channels to use if there are more than one
                excludeEvents: don't include events in the list if they have these IDs
        """
        evts = self[:, self._eventChannels[0]]
        deltaEvt = evts[:-1] - evts[1:]
        deltaIndices = np.where(deltaEvt!=0)[0] #these inds are the LAST entry of old val
        if onsets:
            deltaIndices += 1
        #now get rid of unwanted events
        for unwanted in excludeEvents:
            currentEvents=evts[deltaIndices]
            deltaIndices = deltaIndices[ currentEvents!=unwanted ]
        #now retrieve the events at those indices
        IDs = evts[deltaIndices]
        return deltaIndices, IDs

class ERPs(object):
    """A class with attributes:
        .dat: the raw voltage data
            shape = (nTimepoints, nChannels, nEvents,)
        .conds: the ID of condition for each trial 
            shape = (nEvents,)
        .t: array of times relative to epoch onset 
            shape = (nTimepoints,)
    """
    def __init__(self, dbpa, tStart=-100, tEnd=500, onsets=True, channels=None):
        """:params:
                dbpa: 
                    the DBPA class that provides the data to be epoched
                tStart: 
                    start time relative to epoch onset (in ms)
                tEnd: 
                    end time relative to epoch onset (in ms)
                onsets: 
                    bool, if False then the offsets will be used for epochs
                channels:
                    If None then the data channels from dbpa._dataChannels will be returned
                    Otherwise a list of channels can be specified. e.g. range(122) for everything
        """
        if channels == None:
            channels = dbpa._dataChannels
        #fetch indices and event IDs
        indices, self.IDs = dbpa.getEventIndices(onsets=onsets)
        #get values around those points  in time
        rate = dbpa.samplRate
        nBefore = rate*tStart/1000
        nAfter = rate*tEnd/1000
        self.t = np.linspace(tStart, tEnd, nAfter-nBefore+1)
        #go through those in a loop
        self.dat = np.zeros((len(self.t), len(channels), len(self.IDs)), 'f')
        for evtN, evtIndex in enumerate(indices):
            self.dat[:,:,evtN] = dbpa[ (evtIndex+nBefore):(evtIndex+nAfter+1), channels]
         
        
class DBPAfiltered(EEGfiltered, DBPA):
    """F32 read-only object with included frequency-filteres"""
    def __init__(self,filename,filter_function, **kw_args):
        DBPA.__init__(self,filename,"r+", **kw_args)
        EEGfiltered.__init__(self, filter_function)
        
class DBPAfiltered(DBPA):
    """DBPA read-only object with included frequency-filters"""
    def __init__(self,filename,filter_function):
        DBPA.__init__(self,filename,"r+")
        self._filter_function = filter_function

    @property
    def ff(self):
        return self._filter_function

    def __getitem__(self,item):
        """Calls method of super-class, filters the return value"""
        return self.ff(DBPA.__getitem__(self,item))

class DBPAFilteredWithCache(DBPA):
    """DBPA read-only object with included frequency-filteres"""
    def __init__(self,filename,btype='lp',fl=None,fh=None,border=2,filter_windowed=False,dirname=None):
        from eegpy.ui.eegpylab import freqfilt_eeg
        fd, tmpfn = tempfile.mkstemp(dir=dirname)
        freqfilt_eeg(filename,tmpfn,btype,fl,fh,border,filter_windowed)
        DBPA.__init__(self,tmpfn)

    def __del__(self):
        self.close()
        try:
            os.unlink(self._fn)
        except Exception,e:
            print "Cannot remove file %s"%self._fn, e

if __name__ == "__main__":
    fn = "/home/thorsten/eyecalibrate.000.bh.dat"
    dbpa = DBPA(fn)
